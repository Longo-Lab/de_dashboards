library(tidyverse)
library(data.table)


lab_dir <- file.path('', 'oak', 'stanford', 'scg', 'lab_flongo')
proj_dir <- file.path(lab_dir, 'Tau-PS19_C31_cortex_snRNAseq')
seurat_dir <- file.path(proj_dir, 'seurat_v3')

nameset <- 'PS19_C31'
output_dir <- file.path('.', 'outputs', nameset)  # relative to Rmd


# Get Ensembl genes info
source(file.path(lab_dir, 'de_dashboards', 'scripts', 'genes_info.R'))
ensembl_genes <- get_genes_info(ensembl_ver = '105')


# Generate results
typeouts <- c('MG', 'ASC')
filedates <- c('20220711', '20220721')

for (i in seq(typeouts)) {
  rmarkdown::render(
    input = 'dashboard.Rmd',
    output_file = file.path(output_dir, str_c(typeouts[[i]], '.html')),
    params = list(
      base_dir = seurat_dir,
      ensembl_genes = ensembl_genes,
      nameset = nameset,
      geno = 'PS19',
      drug = 'C31',
      filedate = filedates[[i]],
      typeout = typeouts[[i]]
    )
  )
}
