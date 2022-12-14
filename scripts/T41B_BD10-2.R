library(tidyverse)
library(data.table)


lab_dir <- file.path('', 'oak', 'stanford', 'scg', 'lab_flongo')
proj_dir <- file.path(lab_dir, 'Harry', 'T41B_mapping')
dashboard_dir <- file.path(proj_dir, 'dashboard')

nameset <- 'T41B_BD10-2'
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
    drug = 'BD10-2',
    filedate = '20220914',
    analyses = c('Tg-VvsWt-V', 'Tg-DvsWt-V', 'Tg-DvsTg-V', 'Wt-DvsWt-V'),
    de_cols = c('GENE', 'shrunkenL2FC', 'pvalue', 'svalue', 'log2FoldChange', 'padj'),
    de_names = c('Log2FC (shrunk)', 'Pval', 'Sval', 'Log2FC', 'Pval (adj)'),
    pval_col = 'pvalue',
    lfc_unshrunken_col = 'log2FoldChange',
    fdr_cutoff = 0.05,
    lfc_cutoff_geno = 1,
    lfc_cutoff_drug = 1.5
  )
)
