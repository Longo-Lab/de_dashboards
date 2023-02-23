library(tidyverse)
library(data.table)


lab_dir <- file.path('', 'oak', 'stanford', 'scg', 'lab_flongo')
proj_dir <- file.path(lab_dir, 'PS19_C31_stim')
dashboard_dir <- file.path(proj_dir, 'Gene_level', 'dashboard')

nameset <- 'PS19_C31_stim'
output_dir <- file.path('.', 'outputs', nameset)  # relative to Rmd


# Get Ensembl genes info
source(file.path(lab_dir, 'de_dashboards', 'scripts', 'genes_info.R'))
ensembl_genes <- get_genes_info()


# PS19_C31_stim
proj1 <- 'PS19_C31_stim'

# Generate results
rmarkdown::render(
  input = 'dashboard.Rmd',
  output_file = file.path(output_dir, str_c(proj1, '.html')),
  params = list(
    base_dir = dashboard_dir,
    ensembl_genes = ensembl_genes,
    nameset = proj1,
    geno = 'PS19',
    drug = 'C31',
    analyses = c('Tg-VvsWt-V', 'Tg-DvsWt-V', 'Tg-DvsTg-V', 'Wt-DvsWt-V'),
    de_cols = c('Gene_id', 'log2FoldChange', 'padj', 'pvalue'),
    de_names = c('Log2FC', 'Pval (adj)', 'Pval'),
    pval_col = 'pvalue',
    lfc_unshrunken_col = 'log2FoldChange',
    fdr_cutoff = 0.05,
    lfc_cutoff_geno = 3.5,
    lfc_cutoff_drug = 3.5
  )
)


# PS19_C31_unstim
proj2 <- 'PS19_C31_unstim'

# Generate results
rmarkdown::render(
  input = 'dashboard.Rmd',
  output_file = file.path(output_dir, str_c(proj2, '.html')),
  params = list(
    base_dir = dashboard_dir,
    ensembl_genes = ensembl_genes,
    nameset = proj2,
    geno = 'PS19',
    drug = 'C31',
    analyses = c('Tg-VvsWt-V', 'Tg-DvsWt-V', 'Tg-DvsTg-V', 'Wt-DvsWt-V'),
    de_cols = c('Gene_id', 'log2FoldChange', 'pvalue', 'padj'),
    de_names = c('Log2FC', 'Pval', 'Pval (adj)'),
    pval_col = 'pvalue',
    lfc_unshrunken_col = 'log2FoldChange',
    fdr_cutoff = 0.05,
    lfc_cutoff_geno = 3.5,
    lfc_cutoff_drug = 3.5
  )
)
