import sys, pandas as pd
if len(sys.argv) < 4:
    raise SystemExit("Usage: prep_tableau.py <genus.csv> <metadata.csv> <out.csv>")
genus_path, meta_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

genus = pd.read_csv(genus_path)
meta  = pd.read_csv(meta_path)

keep = [c for c in meta.columns if c.lower() in {
    "sample_id","subject_id","gestational_age","pregnancy_outcome",
    "body_site","visit_number","age","bmi"
}]
meta_slim = meta[keep].copy() if keep else meta

genus_long = genus.melt(id_vars=["Genus"], var_name="sample_id", value_name="abundance")
tab = genus_long.merge(meta_slim, on="sample_id", how="left")
tab.to_csv(out_path, index=False)
print("Exported:", out_path)
