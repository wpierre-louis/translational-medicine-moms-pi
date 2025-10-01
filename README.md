# translational-medicine-moms-pi
Origin
This project uses the Human Microbiome Project 2 (HMP2) MomSPI dataset, a large-scale study characterizing microbial dynamics in pregnancy. The dataset contains multi-omic and clinical metadata across thousands of samples, providing insight into how microbial communities shift in maternal health. 

Focus
While the broader MomSPI dataset includes oral, gut, cervical, and other body sites, this repository focuses on the vaginal microbiome subset during pregnancy. This site-specific lens highlights microbial taxa that play a central role in maternal and neonatal outcomes.

Current Outputs:
- Processed CSV exports filtered to vaginal samples only.

- Data aggregated at multiple levels, mainly: TU, taxonomy, and genus.

- Tableau-ready outputs for straightforward visualization and dashboard building.

Future Outputs:

- Machine learning pipelines (Python) for predictive modeling of microbial community shifts.

- Interactive Jupyter notebooks that demonstrate clustering, ordination, and feature importance analyses.

- Integration of R-based microbiome ecology tools with Python-based ML workflows.

User Guide:

1. Reproduce the vaginal subset exports
```bash 
Rscript scripts/export_momspi_vaginal.R data/processed
```
2. Aggregate taxonomy and prepare Tableau inputs
```bash
Rscript scripts/aggregate_taxa.R \
  data/processed/vaginal_otu.csv \
  data/processed/vaginal_taxonomy.csv \
  data/processed/vaginal_genus.csv

python scripts/prep_tableau.py \
  data/processed/vaginal_genus.csv \
  data/processed/vaginal_metadata.csv \
  data/processed/tableau_vaginal_genus.csv
```

3. Open in Tableau
Load "data/processed/tableau_vaginal_genus.csv" to explore vaginal community composition visually.