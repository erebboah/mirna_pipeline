import sys
import os
import gzip
import regex
import time

# === CONFIGURATION ===

RT_PRIMER = "ACGGGCTAATATTTATCGGTGGAGC"
MAX_MISMATCHES = 6

# 5' barcodes
barcodes_5p = {
    "CAGTCG": "1", "TGACTC": "1", "GCTAGA": "1", "ATCGAT": "1",
    "TCGCAG": "2", "CTCTGA": "2", "AGAGCT": "2", "GATATC": "2",
    "CACGTG": "3", "TGTACC": "3", "GCGTAA": "3", "ATACGT": "3",
    "CAGCGT": "4", "TGCTAC": "4", "GCAGTA": "4", "ATTACG": "4",
    "AGCGTC": "5", "GATCCT": "5", "CTGAAG": "5", "TCATGA": "5",
    "GCATCG": "6", "ATGCTC": "6", "TGCAGA": "6", "CATGAT": "6",
    "GTCGCA": "7", "ACTCTG": "7", "TAGAGC": "7", "CGATAT": "7",
}

# 3' barcodes (after RT primer)
barcodes_3p = {
    "ATCACG": "1", "CGATGT": "2", "TTAGGC": "3", "TGACCA": "4",
    "ACAGTG": "5", "GCCAAT": "6", "CAGATC": "7", "ACTTGA": "8",
    "GATCAG": "9", "TAGCTT": "10", "GGCTAC": "11", "CTTGTA": "12",
    "AGTTCC": "14", "ATGTCA": "15", "GTCCGC": "18", "GTTTCG": "21",
    "CGTACG": "22", "GAGTGG": "23", "GGTAGC": "24", "ACTGAT": "25"
}

# === FUNCTIONS ===

def find_5p_bin(seq):
    prefix = seq[:6]
    for bc, val in barcodes_5p.items():
        if sum(a == b for a, b in zip(prefix, bc)) >= 5:
            return val
    return None

def find_3p_bin(seq):
    match = regex.search(f"({RT_PRIMER}){{e<={MAX_MISMATCHES}}}", seq)
    if match:
        rt_end = match.end()
        bc3_candidate = seq[rt_end:rt_end+6]
        for bc, val in barcodes_3p.items():
            if sum(a == b for a, b in zip(bc3_candidate, bc)) >= 5:
                return val
    return None

def main(fq_path, outdir):
    os.makedirs(outdir, exist_ok=True)
    handles = {}
    read_count = 0
    matched_count = 0

    with gzip.open(fq_path, "rt") as fq:
        while True:
            id_line = fq.readline()
            if not id_line:
                break
            seq = fq.readline().strip()
            plus = fq.readline()
            qual = fq.readline().strip()

            read_count += 1

            bin5 = find_5p_bin(seq)
            if not bin5:
                continue
            bin3 = find_3p_bin(seq)
            if not bin3:
                continue

            matched_count += 1
            pair = f"{bin5}_{bin3}"
            out_path = os.path.join(outdir, f"{pair}.fastq")

            if out_path not in handles:
                handles[out_path] = open(out_path, "w")

            handles[out_path].write(f"{id_line}{seq}\n{plus}{qual}\n")

            if read_count % 100000 == 0:
                print(f"[{time.ctime()}] Processed {read_count:,} reads, matched {matched_count:,}", flush=True)

    for fh in handles.values():
        fh.close()

    print(f"[{time.ctime()}] DONE: Total reads = {read_count:,}, matched = {matched_count:,}", flush=True)

# === ENTRY POINT ===

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python dual_demux.py <input.fastq.gz> <output_dir>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])

