#!/usr/bin/env python3
import argparse, bisect, csv, gzip, itertools, sys
from collections import defaultdict


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_ped(p):
    ped = {}
    for line in opener(p):
        f = line.split()
        if len(f) >= 4 and not line.startswith("#"):
            ped[f[1]] = tuple(x if x not in ("0", "NA", ".") else None for x in f[2:4])
    return ped


def kinship(ped):
    order, seen = [], set()
    def visit(i):
        if i in seen or i is None or i not in ped:
            return
        seen.add(i)
        for p in ped[i]:
            visit(p)
        order.append(i)
    for i in ped:
        visit(i)
    rank = {i: k for k, i in enumerate(order)}
    memo = {}
    def phi(a, b):
        if a is None or b is None:
            return 0.0
        if rank[a] < rank[b]:
            a, b = b, a
        if (a, b) in memo:
            return memo[(a, b)]
        f, m = ped[a]
        v = 0.5 * (1 + phi(f, m)) if a == b else 0.5 * (phi(f, b) + phi(m, b))
        memo[(a, b)] = v
        return v
    return phi


def read_map(p, chrom):
    bp, cm = [], []
    for line in opener(p):
        f = line.split()
        if len(f) >= 4 and f[0] == chrom:
            cm.append(float(f[2]))
            bp.append(int(f[3]))
    def to_cm(x):
        i = bisect.bisect_left(bp, x)
        if i == 0:
            return cm[0]
        if i == len(bp):
            return cm[-1]
        return cm[i - 1] + (cm[i] - cm[i - 1]) * (x - bp[i - 1]) / max(1, bp[i] - bp[i - 1])
    return to_cm


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def intersect(a, b):
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def length(iv, to_cm):
    return sum(to_cm(e) - to_cm(s) for s, e in iv)


def read_iv(p, chrom):
    blocks = []
    with opener(p) as fh:
        for r in csv.DictReader(fh):
            if r["CHROM"] == chrom:
                blocks.append((int(r["start"]), int(r["end"]), {k: v for k, v in r.items() if k not in ("CHROM", "start", "end")}))
    return sorted(blocks, key=lambda b: b[0])


def truth_sets(blocks, ha, hb, mp):
    def label(codes, h):
        c = codes.get(h[0], "")
        if len(c) != 2:
            return None
        return c[(h[1] - 1) ^ (1 if mp else 0)]
    T, N = [], []
    prev = None
    for s, e, codes in blocks:
        la, lb = label(codes, ha), label(codes, hb)
        if la is None or lb is None:
            prev = None
            continue
        st = la == lb
        if prev is not None and prev[0] == st:
            s = prev[1]
        (T if st else N).append((s, e))
        prev = (st, e)
    return merge(T), merge(N)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ibd")
    ap.add_argument("--ped", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--iv")
    ap.add_argument("--mp", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()

    ped = read_ped(a.ped)
    phi = kinship(ped)
    to_cm = read_map(a.map, a.chrom)
    segs = defaultdict(list)
    samples = set()
    with opener(a.ibd) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["chrom"] != a.chrom:
                continue
            ha, hb = (r["sample_a"], int(r["hap_a"])), (r["sample_b"], int(r["hap_b"]))
            samples.update((ha[0], hb[0]))
            segs[tuple(sorted((ha, hb)))].append((int(r["start_bp"]), int(r["end_bp"])))
    segs = {k: merge(v) for k, v in segs.items()}
    blocks = read_iv(a.iv, a.chrom) if a.iv else []
    ivs = set(blocks[0][2]) if blocks else set()
    if blocks:
        lo, hi = blocks[0][0], blocks[-1][1]
    L = to_cm(10 ** 12) - to_cm(0)
    inds = sorted(s for s in samples | (ivs & set(ped)) if s in ped)

    rows = []
    tp = fp = fn = 0.0
    for x, y in itertools.combinations(inds, 2):
        haps = [((x, i), (y, j)) for i in (1, 2) for j in (1, 2)]
        cov = [segs.get(tuple(sorted(h)), []) for h in haps]
        k_obs = sum(length(c, to_cm) for c in cov) / (4 * L)
        ibd2 = merge(intersect(cov[0], cov[3]) + intersect(cov[1], cov[2]))
        row = {"sample_a": x, "sample_b": y, "kinship_exp": round(phi(x, y), 5),
               "kinship_obs": round(k_obs, 5), "ibd1_cM": round(length(merge(sum(cov, [])), to_cm), 3),
               "ibd2_cM": round(length(ibd2, to_cm), 3)}
        if x in ivs and y in ivs:
            t = n = c_tp = c_fp = 0.0
            for h, c in zip(haps, cov):
                T, N = truth_sets(blocks, h[0], h[1], a.mp)
                c = intersect(c, [(lo, hi)])
                t += length(T, to_cm)
                c_tp += length(intersect(c, T), to_cm)
                c_fp += length(intersect(c, N), to_cm)
            row.update(truth_cM=round(t, 3), tp_cM=round(c_tp, 3), fp_cM=round(c_fp, 3),
                       recall=round(c_tp / t, 4) if t else "", precision=round(c_tp / (c_tp + c_fp), 4) if c_tp + c_fp else "")
            tp, fp, fn = tp + c_tp, fp + c_fp, fn + t - c_tp
        rows.append(row)

    keys = list(dict.fromkeys(k for r in rows for k in r))
    out = open(a.out, "w") if a.out else sys.stdout
    w = csv.DictWriter(out, keys, delimiter="\t", restval="")
    w.writeheader()
    w.writerows(rows)
    if rows:
        err = [abs(r["kinship_obs"] - r["kinship_exp"]) for r in rows]
        print(f"pairs={len(rows)} chrom_cM={L:.2f} kinship_MAE={sum(err) / len(err):.4f}", file=sys.stderr)
    if tp + fp + fn:
        print(f"segment_recall={tp / (tp + fn):.4f} segment_precision={tp / (tp + fp) if tp + fp else float('nan'):.4f} "
              f"tp_cM={tp:.2f} fp_cM={fp:.2f} fn_cM={fn:.2f}", file=sys.stderr)


if __name__ == "__main__":
    main()
