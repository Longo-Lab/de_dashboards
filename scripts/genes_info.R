library(data.table)
library(stringr)
library(biomaRt)


get_genes_info <- function(ensembl_set = 'mmusculus_gene_ensembl', ensembl_ver) {
  # Fetch Ensembl genes
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
  lab_dir <- file.path('', 'oak', 'stanford', 'scg', 'lab_flongo')
  ref_dir <- file.path(lab_dir, 'reference', 'biodomains')
  
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
  
  ensembl_genes
}
