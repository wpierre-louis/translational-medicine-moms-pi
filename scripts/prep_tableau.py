#!/usr/bin/env python3
import pandas as pd, sys

if len(sys.argv) < 4:
    print("Usage: python scripts/prep_tableau.py <genus_csv> <metadata_csv> <out_csv> [--join-col <colname>]")
    sys.exit(1)

genus_path, meta_path, out_path = sys.argv[1:4]
join_col = None
if len(sys.argv) > 4 and sys.argv[4] == "--join-col":
    join_col = sys.argv[5]

g = pd.read_csv(genus_path)
m = pd.read_csv(meta_path)
g.columns = [c.strip().strip('"').lower() for c in g.columns]
m.columns = [c.strip().strip('"').lower() for c in m.columns]

def clean(s): return s.astype(str).str.replace("\ufeff","", regex=False).str.strip().str.strip('"').str.strip("'")
m["file_name"] = clean(m["file_name"])

if join_col is None:
    # pick best overlap vs m['file_name']
    candidates = [c for c in ["sample_id","file_name","sample","samplename","id"] if c in g.columns]
    if not candidates:
        candidates = [c for c in g.columns if g[c].dtype == "object"]
    best, best_ov = None, -1
    mset = set(m["file_name"])
    for c in candidates:
        ov = len(set(clean(g[c])) & mset)
        if ov > best_ov:
            best, best_ov = c, ov
    if best is None or best_ov == 0:
        raise SystemExit(f"ERROR: no overlap auto-detected. Pass --join-col <colname>. Genus columns: {list(g.columns)}")
    join_col = best

g[join_col] = clean(g[join_col])
merged = g.merge(m[["file_name","subject_id","visit_number"]],
                 left_on=join_col, right_on="file_name", how="left")

if join_col != "file_name":
    merged = merged.rename(columns={join_col: "sample_id"})
else:
    merged = merged.rename(columns={"file_name": "sample_id"})
if "file_name" in merged.columns and "sample_id" in merged.columns and merged["file_name"].equals(merged["sample_id"]):
    merged = merged.drop(columns=["file_name"])

merged.to_csv(out_path, index=False)
print(f"Used join column: {join_col} | Exported: {out_path}")
