library(tidyverse)
library(data.table)


lab_dir <- file.path('', 'oak', 'stanford', 'scg', 'lab_flongo')
proj_dir <- file.path(lab_dir, 't41b_BD10-2_stim')
dashboard_dir <- file.path(proj_dir, 'dashboard')

nameset <- 'T41B_BD10-2_stim'
output_dir <- file.path('.', 'outputs', nameset)  # relative to Rmd


# Get Ensembl genes info
source(file.path(lab_dir, 'de_dashboards', 'scripts', 'genes_info.R'))
ensembl_genes <- get_genes_info()


# T41B_BD10-2_stim & T41B_BD10-2_unstim (same settings)
projs <- c('T41B_BD10-2_stim', 'T41B_BD10-2_unstim')

# Generate results
for (i in projs) {
  rmarkdown::render(
    input = 'dashboard.Rmd',
    output_file = file.path(output_dir, str_c(i, '.html')),
    params = list(
      base_dir = dashboard_dir,
      ensembl_genes = ensembl_genes,
      nameset = i,
      geno = 'APPL/S',
      drug = 'BD10-2',
      de_cols = c('Gene_id', 'log2FoldChange', 'pvalue', 'padj'),
      de_names = c('Log2FC', 'Pval', 'Pval (adj)'),
      pval_col = 'pvalue',
      lfc_unshrunken_col = 'log2FoldChange',
      fdr_cutoff = 0.05,
      lfc_cutoff_geno = 2.5,
      lfc_cutoff_drug = 1
    )
  )
}


# T41B_stim
proj3 <- 'T41B_stim'

# Generate results
rmarkdown::render(
  input = 'dashboard.Rmd',
  output_file = file.path(output_dir, str_c(proj3, '.html')),
  params = list(
    base_dir = dashboard_dir,
    ensembl_genes = ensembl_genes,
    nameset = proj3,
    geno = 'APPL/S',
    drug = 'TBS',
    analyses = c('Tg-VvsWt-V', 'Tg-DvsWt-V', 'Tg-DvsTg-V', 'Wt-DvsWt-V'),
    de_cols = c('Gene_id', 'log2FoldChange', 'padj', 'pvalue'),
    de_names = c('Log2FC', 'Pval (adj)', 'Pval'),
    pval_col = 'pvalue',
    lfc_unshrunken_col = 'log2FoldChange',
    fdr_cutoff = 0.05,
    lfc_cutoff_geno = 3,
    lfc_cutoff_drug = 5
  )
)
