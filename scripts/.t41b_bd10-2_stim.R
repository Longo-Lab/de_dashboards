library(tidyverse)
library(data.table)
library(biomaRt)


lab_dir <- file.path('', 'oak', 'stanford', 'scg', 'lab_flongo')
proj_dir <- file.path(lab_dir, 't41b_BD10-2_stim')
seurat_dir <- file.path(proj_dir, 'dashboard')
ref_dir <- file.path(lab_dir, 'reference', 'biodomains')

output_dir <- file.path('.', 'outputs', 't41b_BD10-2_stim')  # relative to Rmd
templates_dir <- file.path('.', 'templates')  # relative to Rmd


# Fetch Ensembl genes
ensembl_set <- 'mmusculus_gene_ensembl'
ensembl_ver <- '105'

ensembl <- useEnsembl('ensembl', dataset = ensembl_set, version = ensembl_ver)
ensembl_genes <- data.table(getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name'),
  filters = 'chromosome_name',
  values = c(1:19, 'X', 'Y', 'MT'),
  mart = ensembl
))

keep_categ <- grep('pseudo', unique(ensembl_genes$gene_biotype), value = T, invert = T)
keep_categ <- grep('TEC', keep_categ, value = T, invert = T)
ensembl_genes <- ensembl_genes[gene_biotype %in% keep_categ]

# Read in modules
modules <- c('TREAT-AD.Mouse.Genes.txt', 'TMT-LP.Mouse.Genes.txt', 'Mostafavi_etal.Mouse.Genes.txt', 'Milind_etal.Mouse.Genes.txt', 'AD_modules_full.tsv')
names(modules) <- c('treatAD', 'tmtAD', 'mostafavi', 'milind', 'wan')

for (m in names(modules)) {
  module <- fread(file.path(ref_dir, modules[[m]]))
  setnames(module, c('Mmus_GeneID', 'Module'), c('GENE', 'module'), skip_absent = T)
  
  # Collapse by GENE (Ensembl ID)
  module <- module[, setNames(.(str_flatten(unique(module), ',')), m), by = GENE]
  
  # Merge modules
  ensembl_genes <- merge(ensembl_genes, module, by.x = 'ensembl_gene_id', by.y = 'GENE', all.x = T)
}


# Generate results
subclasses <- c('t41b_BD10-2_stim')
filedates <- c('20220911')

for (i in seq(subclasses)) {
  rmarkdown::render(
    input = 'bulk_dashboard.Rmd',
    output_file = file.path(output_dir, str_c('t41b_BD10-2_stim', '.html')),
    params = list(
      subclass = subclasses[[i]],
      filedate = filedates[[i]],
      seurat_dir = seurat_dir,
      templates_dir = templates_dir,
      ensembl_genes = ensembl_genes,
      nameset = 't41b_BD10-2_stim',
      geno = 'APPL/S',
      drug = 'BD10-2'
    )
  )
}
