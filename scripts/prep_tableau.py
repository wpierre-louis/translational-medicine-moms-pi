import sys
import pandas as pd

if len(sys.argv) != 4:
    print(__doc__.strip()); sys.exit(1)

genus_path, meta_path, out_path = sys.argv[1:4]

# -------- Load --------
g = pd.read_csv(genus_path)
m = pd.read_csv(meta_path)

# -------- Helpers --------
def ci_lookup(cols, target):
    """Return the actual column name in `cols` that matches `target` case-insensitively."""
    tl = target.lower()
    for c in cols:
        if c.lower() == tl:
            return c
    return None

def clean(s: pd.Series) -> pd.Series:
    return (
        s.astype(str)
         .str.replace("\ufeff", "", regex=False)
         .str.strip().str.strip('"').str.strip("'")
    )

# -------- Identify required columns (case-insensitive) --------
genus_col = ci_lookup(g.columns, "Genus")
if genus_col is None:
    raise SystemExit('ERROR: Genus column not found in genus CSV (expected "Genus").')

file_name_col = ci_lookup(m.columns, "file_name")
if file_name_col is None:
    raise SystemExit('ERROR: "file_name" column not found in metadata CSV.')

# Optional metadata fields
opt_fields = ["subject_id", "visit_number", "subject_gender", "subject_race"]
meta_keep = [c for name in [file_name_col] + opt_fields for c in m.columns if c.lower() == name.lower()]

# Clean the metadata join key
m[file_name_col] = clean(m[file_name_col])

# -------- Wide → Long (keep *exact* "Genus" column name) --------
value_cols = [c for c in g.columns if c != genus_col]
if not value_cols:
    raise SystemExit("ERROR: genus CSV appears to have no sample columns to melt.")

g_long = g.melt(id_vars=[genus_col], var_name="sample_id", value_name="count")

# Coerce counts to numeric
g_long["count"] = pd.to_numeric(g_long["count"], errors="coerce").fillna(0)

# -------- Merge metadata --------
merged = (
    g_long.merge(m[meta_keep], left_on="sample_id", right_on=file_name_col, how="left")
          .drop(columns=[file_name_col])
)
# -------- Write --------
merged.to_csv(out_path, index=False)

# -------- Log --------
n_samples = merged["sample_id"].nunique()
n_features = merged[genus_col].nunique()
print(f"Exported: {out_path} | rows={len(merged)} | samples={n_samples} | genera={n_features}")
