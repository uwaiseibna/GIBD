#!/usr/bin/env python3
import csv, os, subprocess, sys, tempfile

here = os.path.dirname(os.path.abspath(__file__))
d = tempfile.mkdtemp()
M = 10 ** 6
open(f"{d}/f.ped", "w").write("F\tP\tNA\tNA\t1\t0\nF\tQ\tNA\tNA\t2\t0\nF\tK1\tP\tQ\t1\t0\nF\tK2\tP\tQ\t2\t0\nF\tG\tK1\tNA\t1\t0\n")
open(f"{d}/m.map", "w").write("".join(f"c1\t.\t{i}\t{i * M}\n" for i in range(101)))
open(f"{d}/iv.csv", "w").write("CHROM,start,end,K1,K2,P,Q\nc1,0,40000000,AC,AD,AB,CD\nc1,41000000,100000000,AC,BD,AB,CD\n")
seg = [("P", 1, "K1", 1, 0, 100), ("Q", 1, "K1", 2, 0, 100), ("P", 1, "K2", 1, 0, 40), ("P", 2, "K2", 1, 41, 100),
       ("Q", 2, "K2", 2, 0, 100), ("K1", 1, "K2", 1, 0, 30), ("K1", 2, "K2", 2, 50, 55), ("G", 1, "K1", 1, 0, 50)]
with open(f"{d}/g.tsv", "w") as f:
    f.write("sample_a\thap_a\tsample_b\thap_b\tchrom\tstart_bp\tend_bp\tcM\ttype\tinverted\tinsertion\tdeletion\tinv_allele\n")
    for a, ha, b, hb, s, e in seg:
        f.write(f"{a}\t{ha}\t{b}\t{hb}\tc1\t{s * M}\t{e * M}\t{e - s}\tIBD\t0\t0\t0\t.\n")
r = subprocess.run([sys.executable, f"{here}/score_pedigree.py", f"{d}/g.tsv", "--ped", f"{d}/f.ped", "--map", f"{d}/m.map",
                    "--chrom", "c1", "--iv", f"{d}/iv.csv", "--out", f"{d}/o.tsv"], capture_output=True, text=True, check=True)
rows = {(x["sample_a"], x["sample_b"]): x for x in csv.DictReader(open(f"{d}/o.tsv"), delimiter="\t")}
k = rows[("K1", "K2")]
checks = [
    (float(k["kinship_exp"]), 0.25), (float(rows[("P", "Q")]["kinship_exp"]), 0.0),
    (float(rows[("K1", "P")]["kinship_exp"]), 0.25), (float(rows[("G", "K1")]["kinship_obs"]), 0.125),
    (float(k["truth_cM"]), 40.0), (float(k["tp_cM"]), 30.0), (float(k["fp_cM"]), 5.0),
    (float(k["recall"]), 0.75), (float(k["precision"]), round(30 / 35, 4)),
    (float(rows[("K2", "P")]["truth_cM"]), 99.0), (float(rows[("K2", "P")]["recall"]), 1.0),
    (float(rows[("K1", "Q")]["kinship_obs"]), 0.25), (float(k["kinship_obs"]), 0.0875),
    (float(rows[("K2", "Q")]["fp_cM"]), 0.0), (float(rows[("G", "P")]["kinship_exp"]), 0.125),
    (float(rows[("G", "K2")]["kinship_exp"]), 0.125), (float(rows[("G", "Q")]["kinship_obs"]), 0.0),
]
bad = [(i, a, b) for i, (a, b) in enumerate(checks) if abs(a - b) > 1e-6]
print("PASS score_pedigree" if not bad else f"FAIL score_pedigree {bad}")
print(r.stderr.strip())
sys.exit(1 if bad else 0)
