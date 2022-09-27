library(tidyverse)
library(data.table)


lab_dir <- file.path('', 'oak', 'stanford', 'scg', 'lab_flongo')
proj_dir <- file.path(lab_dir, 'Harry', 'T41B_mapping')
dashboard_dir <- file.path(proj_dir, 'dashboard')

nameset <- 'T41B_BD10-2'
output_dir <- file.path('.', 'outputs', nameset)  # relative to Rmd


# Get Ensembl genes info
source(file.path(lab_dir, 'de_dashboards', 'scripts', 'genes_info.R'))
ensembl_genes <- get_genes_info(ensembl_ver = '105')


# Generate results
rmarkdown::render(
  input = 'dashboard.Rmd',
  output_file = file.path(output_dir, str_c(nameset, '.html')),
  params = list(
    base_dir = dashboard_dir,
    ensembl_genes = ensembl_genes,
    nameset = nameset,
    geno = 'APPL/S',
    drug = 'BD10-2',
    filedate = '20220914',
    de_cols = c('GENE', 'shrunkenL2FC', 'svalue', 'log2FoldChange', 'padj', 'pvalue'),
    de_names = c('Log2FC (shrunk)', 'Sval', 'Log2FC', 'Pval (adj)', 'Pval'),
    pval_col = 'pvalue',
    lfc_unshrunken_col = 'log2FoldChange',
    fdr_cutoff = 0.5,
    lfc_cutoff_geno = 1,
    lfc_cutoff_drug = 1.5
  )
)