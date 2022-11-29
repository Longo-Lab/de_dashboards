library(tidyverse)
library(data.table)


lab_dir <- file.path('', 'oak', 'stanford', 'scg', 'lab_flongo')
proj_dir <- file.path(lab_dir, 't41b_BD10-2_stim')
dashboard_dir <- file.path(proj_dir, 'dashboard')

nameset <- 'T41B_stim'
output_dir <- file.path('.', 'outputs', nameset)  # relative to Rmd


# Get Ensembl genes info
source(file.path(lab_dir, 'de_dashboards', 'scripts', 'genes_info.R'))
ensembl_genes <- get_genes_info()


# Generate results
rmarkdown::render(
  input = 'dashboard.Rmd',
  output_file = file.path(output_dir, str_c(nameset, '.html')),
  params = list(
    base_dir = dashboard_dir,
    ensembl_genes = ensembl_genes,
    nameset = nameset,
    geno = 'APPL/S',
    drug = 'TBS',
    filedate = '20220922',
    de_cols = c('Gene_id', 'log2FoldChange', 'padj', 'pvalue'),
    de_names = c('Log2FC', 'Pval (adj)', 'Pval'),
    pval_col = 'pvalue',
    lfc_unshrunken_col = 'log2FoldChange',
    fdr_cutoff = 0.05,
    lfc_cutoff_geno = 3,
    lfc_cutoff_drug = 5
  )
)
