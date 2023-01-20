library(tidyverse)
library(data.table)


lab_dir <- file.path('', 'oak', 'stanford', 'scg', 'lab_flongo')
proj_dir <- file.path(lab_dir, 'Tau-PS19_C31_cortex_snRNAseq')
seurat_dir <- file.path(proj_dir, 'seurat_v3')

nameset <- 'PS19_C31'
output_dir <- file.path('.', 'outputs', nameset)  # relative to Rmd


# Get Ensembl genes info
source(file.path(lab_dir, 'de_dashboards', 'scripts', 'genes_info.R'))
ensembl_genes <- get_genes_info()


# Generate results
typeouts <- c('MG', 'ASC', 'OPC', 'L4_IT_CTX')

for (i in typeouts) {
  rmarkdown::render(
    input = 'dashboard.Rmd',
    output_file = file.path(output_dir, str_c(i, '.html')),
    params = list(
      base_dir = seurat_dir,
      ensembl_genes = ensembl_genes,
      nameset = nameset,
      geno = 'PS19',
      drug = 'C31',
      round_num = 'R6',
      typeout = i,
      analyses = c('Tg-VvsWt-V', 'Tg-DvsWt-V', 'Tg-DvsTg-V', 'Wt-DvsWt-V')
    )
  )
}
