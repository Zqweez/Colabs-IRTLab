"""
Quality-aware Sanger pipeline (AB1 → trimmed FASTA → Reverse complement → MAFFT → consensus → combined consensus).
Uses MAFFT in subprocess conda install -c bioconda mafft
"""

import os, re, subprocess
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from colorama import Fore, Style

# Configuration
print(Fore.CYAN + "\nStarting Sanger sequencing pipeline at", datetime.now().strftime("%Y-%m-%d %H:%M:%S") + Style.RESET_ALL)
# Directory with raw chromatograms
ab1_dir = Path("data/sanger/00_ab1_raw_reads")
    
# Output folders
trimmed_dir    = Path("data/sanger/01_ab1_trimmed_reads")
aligned_dir    = Path("data/sanger/02_ab1_aligned_read_pairs")
consensus_dir  = Path("data/sanger/03_ab1_assembled_sequences")

for p in (trimmed_dir, aligned_dir, consensus_dir):
    p.mkdir(parents=True, exist_ok=True)

# Regex for pairing
name_re = re.compile(r"^(?P<base>.+?)[_-]?(?P<orient>[FR])(?:-?PREMIX)?\.ab1$", re.IGNORECASE)

# Quality-trimming parameters
WINDOW      = 20   # sliding-window size
CUTOFF      = 20   # require mean Q >= CUTOFF inside the window
MIN_LEN     = 200   # discard reads trimmed shorter than this
MISMATCH_N  = 3    # if |Q1-Q2| < MISMATCH_N and bases differ -> call N

# Find paired ab1 files
sample_map = defaultdict(dict)   # {sample : {'F': Path, 'R': Path}}

for f in ab1_dir.iterdir():
    m = name_re.match(f.name)
    if not m:
        continue
    base   = m.group('base')
    orient = m.group('orient').upper()
    sample_map[base][orient] = f

pairs = [
    {"name": s, "fwd": v['F'], "rev": v['R']}
    for s, v in sample_map.items() if {'F', 'R'} <= v.keys()
]

if not pairs:
    raise SystemExit("No complete F/R pairs found.")

print(f"Found {len(pairs)} chromatogram pairs")

# Quality-aware trimming function
def trim_by_quality(name, seq_record, window=WINDOW, cutoff=CUTOFF):
    """
    Returns (trimmed_seq_record, trimmed_qualities) or (None, None) if too short.
    """
    quals = seq_record.letter_annotations["phred_quality"]
    seq   = str(seq_record.seq)

    # 5' trim
    start = 0
    for i in range(len(seq) - window + 1):
        if sum(quals[i:i+window]) / window >= cutoff:
            start = i
            break

    # 3' trim
    end = len(seq)
    for i in range(len(seq) - window, -1, -1):
        if sum(quals[i:i+window]) / window >= cutoff:
            end = i + window
            break

    trimmed_seq   = seq[start:end]
    trimmed_qual  = quals[start:end]

    if len(trimmed_seq) < MIN_LEN:
        return None, None

    tr_record = SeqRecord(Seq(trimmed_seq),
                          id=name,
                          description=f"ID = {seq_record.id} trimmed Q>{cutoff}")
    tr_record.letter_annotations["phred_quality"] = trimmed_qual
    return tr_record, trimmed_qual

# Consensus generation from alignment
def consensus_from_alignment(aln, qual1, qual2, sample):
    """
    aln   : Biopython MultipleSeqAlignment with 2 seqs (F, R-revcomp)
    qual1 : list of ints for aligned[0] (gaps ignored by index mapping)
    qual2 : list of ints for aligned[1]
    Returns SeqRecord consensus.
    """
    if len(aln) != 2:
        raise ValueError(f"{sample}: alignment does not have 2 sequences")

    s1, s2 = str(aln[0].seq), str(aln[1].seq)

    idx1 = idx2 = 0
    out  = []

    for b1, b2 in zip(s1, s2):
        q1 = qual1[idx1] if b1 != '-' else -1
        q2 = qual2[idx2] if b2 != '-' else -1

        if b1 != '-': idx1 += 1
        if b2 != '-': idx2 += 1

        if b1 == '-':
            out.append(b2)
        elif b2 == '-':
            out.append(b1)
        elif b1 == b2:
            out.append(b1)
        else:
            # choose higher quality; if similar, call N
            if abs(q1 - q2) < MISMATCH_N:
                out.append('N')
            else:
                out.append(b1 if q1 > q2 else b2)

    seq = ''.join(out).strip('N-')
    return SeqRecord(Seq(seq), id=sample, description=f"consensus len={len(seq)}")

# Main processing loop
all_consensus = []

for p in pairs:
    name = p['name']
    print(f"\n▶ Processing {name}")

    # --- read chromatograms ---
    f_rec = SeqIO.read(p['fwd'], "abi")
    r_rec = SeqIO.read(p['rev'], "abi")

    # --- trim reads ---
    f_trim, fq = trim_by_quality(name, f_rec)
    r_trim, rq = trim_by_quality(name, r_rec)

    # ───────────── decision tree ──────────────────────────────────────
    if not f_trim and not r_trim:
        print(Fore.RED + "  ✖  Both reads failed QC; skipping." + Style.RESET_ALL)
        continue              # nothing to do

    elif f_trim and r_trim:
        # reverse-complement the reverse read AFTER trimming
        r_trim_rc = SeqRecord(r_trim.seq.reverse_complement(),
                                id=r_trim.id + "_revcomp",
                                description="revcomp")
        rq = rq[::-1]  # reverse qualities to match rev-comp orientation

        # --- write trimmed FASTA for record-keeping ---
        f_trim_fn = trimmed_dir / f"{name}_F.trim.fasta"
        r_trim_fn = trimmed_dir / f"{name}_R.trim.fasta"
        SeqIO.write(f_trim, f_trim_fn, "fasta")
        SeqIO.write(r_trim_rc, r_trim_fn, "fasta")

        # --- align with MAFFT ---
        temp_in  = aligned_dir / f"{name}_in.fasta"
        aln_out  = aligned_dir / f"{name}_aligned.fasta"
        SeqIO.write([f_trim, r_trim_rc], temp_in, "fasta")

        mafft_cmd = ["mafft", "--auto", "--quiet", str(temp_in)]
        with open(aln_out, "w") as out_f:
            subprocess.run(mafft_cmd, stdout=out_f, check=True)

        aln = AlignIO.read(aln_out, "fasta")
        print(Fore.GREEN + f"  ✓ Aligned with MAFFT" + Style.RESET_ALL)

        # --- consensus ---
        cons = consensus_from_alignment(aln, fq, rq, name)
        cons_fn = consensus_dir / f"{name}_consensus.fasta"
        SeqIO.write(cons, cons_fn, "fasta")
        all_consensus.append(cons)
        print(Fore.CYAN + f"  ↳ Consensus saved → {cons_fn.name} (len {len(cons.seq)})" + Style.RESET_ALL)

    else:
        # Only one read passed QC, use it as is
        good_rec  = f_trim or r_trim           # whichever is not None
        quals     = fq if f_trim else rq
        # if reverse read was the keeper, reverse-complement it
        if good_rec is r_trim:
            good_rec.seq  = good_rec.seq.reverse_complement()
            good_rec.id  += "_revcomp"
            quals = quals[::-1]

        # save trimmed FASTA
        single_fn = trimmed_dir / f"{name}_single.trim.fasta"
        SeqIO.write(good_rec, single_fn, "fasta")

        # single-read “consensus” = the read itself
        cons      = SeqRecord(good_rec.seq,
                            id=name,
                            description=f"single-read consensus len={len(good_rec)}")
        cons_fn   = consensus_dir / f"{name}_consensus.fasta"
        SeqIO.write(cons, cons_fn, "fasta")
        all_consensus.append(cons)
        print(Fore.YELLOW + f"  ↳ Single-read consensus saved → {cons_fn.name}" + Style.RESET_ALL)

# Make a combined FASTA of all consensus sequences
if all_consensus:
    master = consensus_dir / "ALL_SAMPLES_consensus.fasta"
    SeqIO.write(all_consensus, master, "fasta")
    print(f"\n Combined FASTA written to {master}")
else:
    print("\nNo consensus sequences generated.")

print(Fore.CYAN + "\n Pipeline complete at", datetime.now().strftime("%Y-%m-%d %H:%M:%S") + Style.RESET_ALL)
