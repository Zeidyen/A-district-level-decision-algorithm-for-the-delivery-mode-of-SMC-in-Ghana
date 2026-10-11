"""Compare a fresh pipeline run (data/derived) with a previous one (data/derived_prev_*).

Usage: python3 scripts/compare_runs.py data/derived_prev_YYYYMMDD_HHMM
Reports, for every CSV in both runs: identical, or the largest numeric change per column
and any change in text columns (e.g. district class). Files only in one run are listed.
"""
import sys, os
import pandas as pd

old, new = sys.argv[1], "data/derived"
fo = {f for f in os.listdir(old) if f.endswith(".csv")}
fn = {f for f in os.listdir(new) if f.endswith(".csv")}
for f in sorted(fo - fn):
    print(f"ONLY IN OLD  {f}")
for f in sorted(fn - fo):
    print(f"NEW FILE     {f}")
for f in sorted(fo & fn):
    a, b = pd.read_csv(os.path.join(old, f)), pd.read_csv(os.path.join(new, f))
    if a.equals(b):
        print(f"identical    {f}")
        continue
    if a.shape != b.shape or list(a.columns) != list(b.columns):
        print(f"CHANGED      {f}: shape {a.shape} -> {b.shape}; columns equal: {list(a.columns) == list(b.columns)}")
        continue
    keys = [c for c in ("region", "Region", "district", "District", "scenario", "k", "cluster", "term", "model", "subset", "rule", "mode", "band", "measure")
            if c in a.columns]
    if keys:
        a, b = a.sort_values(keys).reset_index(drop=True), b.sort_values(keys).reset_index(drop=True)
    out = []
    for c in a.columns:
        if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]):
            d = (a[c] - b[c]).abs().max()
            if d > 0:
                out.append(f"{c} max|diff| {d:.4g}")
        else:
            n = int((a[c].astype(str) != b[c].astype(str)).sum())
            if n:
                out.append(f"{c}: {n} rows differ")
    print(f"CHANGED      {f}: " + ("; ".join(out) if out else "row order only"))
