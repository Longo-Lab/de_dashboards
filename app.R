library(tidyverse)
library(data.table)
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(gprofiler2)


# ----------------
# Define settings
# ----------------

lab_dir <- file.path('', 'oak', 'stanford', 'scg', 'lab_flongo')
proj_dir <- file.path(lab_dir, 'Tau-PS19_C31_cortex_snRNAseq')
seurat_dir <- file.path(proj_dir, 'seurat_v3')

base_dir <- seurat_dir
nameset <- 'PS19_C31'
geno <- 'PS19'
drug <- 'C31'
round_num <- 'R6'
analyses <- c('Tg-VvsWt-V', 'Tg-DvsWt-V', 'Tg-DvsTg-V', 'Wt-DvsWt-V')
de_cols <- c('GENE', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2')
de_names <- c('Log2FC', 'Pval (adj)', 'Pct 1', 'Pct 2')
pval_col <- 'p_val'
lfc_unshrunken_col <- 'avg_log2FC'
fdr_cutoff <- 0.05
lfc_cutoff_geno <- 0.55 
lfc_cutoff_drug <- 0.5
is_sc <- T
term_size_limit <- 2000

# Add project directories to be able to embed images
img_dir <- 'imgs'
addResourcePath(img_dir, seurat_dir)

# Cell type selection
cell_types <- c(
  'MG' = 'MG',
  'ASC' = 'ASC',
  'L4/5 IT CTX' = 'L4_IT_CTX'
)


# ------------------------
# Define global variables
# ------------------------

# Load genes data
ensembl_genes <- readRDS('www/ensembl_genes.RDS')

# Biodomain modules
modules <- c('TREAT-AD.Mouse.Genes.txt', 'TMT-LP.Mouse.Genes.txt', 'Mostafavi_etal.Mouse.Genes.txt', 'Milind_etal.Mouse.Genes.txt', 'AD_modules_full.tsv')
names(modules) <- c('treatAD', 'tmtAD', 'mostafavi', 'milind', 'wan')

# WT columns
include_wt <- ifelse(length(analyses) == 4, T, F)

# Gene ID variable
ensembl_var <- ifelse(is_sc, 'external_gene_name', 'ensembl_gene_id')

# Genes columns
gene_col <- de_cols[[1]]
lfc_col <- de_cols[[2]]
p_col <- de_cols[[3]]

lfc_geno <- str_c(lfc_col, 'geno', sep = '.')
lfc_drug <- str_c(lfc_col, 'drug', sep = '.')
lfc_geno_drug <- str_c(lfc_col, 'geno_drug', sep = '.')
p_geno <- str_c(p_col, 'geno', sep = '.')
p_drug <- str_c(p_col, 'drug', sep = '.')
p_geno_drug <- str_c(p_col, 'geno_drug', sep = '.')

# Comparisons columns
var_names <- c('geno', 'geno_drug', 'drug')

if (include_wt) {
  var_names <- c(var_names, 'drug_wt')
}

names(analyses) <- var_names

# Genes list table
analyses_cols <- list(
  htmltools::tags$th(colspan = length(de_names), style = 'background-color:#f3f7eb;border-bottom:none;text-align: center;', geno),
  htmltools::tags$th(colspan = length(de_names), style = 'background-color:#fdf2f1;border-bottom:none;text-align: center;', str_c(geno, drug, sep = '_')),
  htmltools::tags$th(colspan = length(de_names), style = 'background-color:#eef8f9;border-bottom:none;text-align: center;', drug)
)

if (include_wt) {
  analyses_cols[[4]] <- htmltools::tags$th(colspan = length(de_names), style = 'background-color:#f9f2ff;border-bottom:none;text-align: center;', str_c(drug, 'wt', sep = '_'))
}

sketch <- htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(rowspan = 2, 'Ensembl ID'),
      th(rowspan = 2, 'Symbol'),
      th(rowspan = 2, 'Category'),
      th(rowspan = 2, 'Direct'),
      analyses_cols,
      th(colspan = 5, style = 'background-color:#eeeeee;border-bottom:none;text-align: center;', 'Modules'),
      th(colspan = 2, style = 'background-color:#f9f3d7;border-bottom:none;text-align: center;', 'TF')
    ),
    tr(
      lapply(rep(de_names, length(analyses_cols)), th),
      th('TreatAD'),
      th('TmtAD'),
      th('Mostafavi'),
      th('Milind'),
      th('Wan'),
      th('AnimalTFDB'),
      th('ChEMBL')
    )
  )
))

modules_col <- 4 + length(analyses) * length(de_names)

render_tooltip <- sprintf("
  function(row, data) {
    // 2nd column (Ensembl gene info)
    var info = data[1].split('|'),
        title = 'Chromosome: ' + info[2] + '&#010;Biotype: ' + info[1];
    
    $('td:eq(1)', row)
      .html('<span class=\"tbl-info\" title=\"' + title + '\">' + info[0] + '</span>');
      
    // Module columns
    var m = %s;
    for (var i = m; i < m + 5; i++) {
      var modules = data[i] ? data[i].split('|')[1].split(',') : [];
      $('td:eq(' + i + ')', row)
        .html('<span class=\"tbl-info\" title=\"' + modules.join('&#010;') + '\">' + modules.length + '</span>');
    }
    
    // TF columns
    $('td:eq(' + (m + 5) + ')', row)
      .html(data[m + 5] === null ? data[m + 5] : data[m + 5].slice(0, -3));
      
    var c = '';
    if (data[m + 6] !== null) {
      var info = data[m + 6].split('~'),
          num = parseInt(info[0]),
          labels = info[1].split('|'),
          source = info[2];
    
      c = '<span class=\"tbl-info\" title=\"Source: ' + source + '&#010;&#010;' + labels.join('&#010;') + '\">' + num + '</span>';
    }
      
    $('td:eq(' + (m + 6) + ')', row)
      .html(c);
  }
  ", modules_col)


# -----------------
# Define functions
# -----------------

get_clusters <- function(typeout) {
  print(str_c('get_clusters() for ', typeout, '...'))
  
  # Obtain clusters
  if (is_sc) {
    file_regex <- str_c(round_num, '\\.', nameset, '\\.(.+)\\.RNA.+\\.csv\\.gz')
    file_matches <- list.files(file.path(base_dir, typeout), pattern = file_regex)
    clusters <- unique(str_match(file_matches, file_regex)[,2])
    clusters <- append(clusters[clusters != typeout], typeout, after = 0)  # order subclass before clusters
  } else {
    clusters <- nameset
  }
  
  names(clusters) <- clusters
  
  clusters
}

get_genes_files <- function(typeout, clusters) {
  print(str_c('get_gene_files() for ', typeout, '...'))
  
  # Read in data
  lapply(clusters, function(cl) lapply(analyses, function(a) {
    fn <- ifelse(is_sc, str_c(nameset, '.', cl), nameset)
    
    fread(file.path(base_dir, typeout, str_c(round_num, fn, 'RNA', a, 'csv.gz', sep = '.')))
  }))
}

get_results <- function(clusters, genes_files) {
  print(str_c('get_results() for ', names(clusters)[[1]], '...'))
  
  # Merge data
  lapply(clusters, function(cl) {
    
    result <- genes_files[[cl]][[1]][!(is.na(get(gene_col))), ..de_cols]
    
    for (i in head(seq(genes_files[[cl]]), -1)) {
      gene_file <- genes_files[[cl]][[i + 1]][, ..de_cols]
      result <- merge(result, gene_file, by = gene_col, all = T, suffixes = c('', str_c('.', names(analyses)[[i + 1]])))
    }
    
    result <- merge(ensembl_genes, result, by.x = ensembl_var, by.y = gene_col)
    setnames(result, de_cols, str_c(de_cols, '.', names(analyses)[[1]]), skip_absent = T)
    
    # Add gene category
    result[, category := fcase(
      (get(lfc_geno) > 0 & get(p_geno) < fdr_cutoff) & ((get(lfc_geno_drug) < 0 & get(p_geno_drug) < fdr_cutoff) | (get(lfc_drug) < 0 & get(p_drug) < fdr_cutoff)), 'suppression',
      (get(lfc_geno) > 0 & get(p_geno) < fdr_cutoff) & (get(lfc_drug) > 0 & get(p_drug) < fdr_cutoff), 'compensatory enhancement',
      (get(lfc_geno) < 0 & get(p_geno) < fdr_cutoff) & ((get(lfc_geno_drug) > 0 & get(p_geno_drug) < fdr_cutoff) | (get(lfc_drug) > 0 & get(p_drug) < fdr_cutoff)), 'enhancement',
      (get(lfc_geno) < 0 & get(p_geno) < fdr_cutoff) & (get(lfc_drug) < 0 & get(p_drug) < fdr_cutoff), 'compensatory suppression',
      default = NA_character_
    )]
    
    # Add gene direction
    result[, direct := fcase(
      category == 'suppression' & (get(lfc_drug) < 0 & get(p_drug) < fdr_cutoff), T,
      category == 'enhancement' & (get(lfc_drug) > 0 & get(p_drug) < fdr_cutoff), T,
      category %in% c('suppression', 'enhancement'), F,
      default = NA
    )]
    
    # Rearrange columns
    var_names <- c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name', 'category', 'direct')
    tf_names <- c('Family', 'is_mouse_chembl', 'chembl_protein_labels')
    setcolorder(
      result, 
      c(
        var_names,
        setdiff(names(result), c(var_names, names(modules), tf_names)),
        names(modules),
        tf_names
      )
    )
  })
}

add_count <- function(x, s1 = ',', s2 = '|') {
  str_c(str_pad(str_count(x, s1) + 1, width = 2, side = 'left', pad = '0'), x, sep = s2)
}

get_gprofiler <- function(gp_e, gp_s, gp_n, gp_a) {
  print(str_c('gProfiler for ', str_trim(str_c(gp_n, ' enhancement/suppression...'))))
  
  e_title <- str_to_sentence(str_trim(str_c(gp_n, ' enhancement')))
  s_title <- str_to_sentence(str_trim(str_c(gp_n, ' suppression')))
  
  e_genes <- gp_e$external_gene_name
  s_genes <- gp_s$external_gene_name
  
  # gProfiler plot
  if (!is_empty(e_genes) && !is_empty(s_genes)) {
    g <- gost(
      query = list(
        e = e_genes,
        s = s_genes
      ),
      organism = 'mmusculus',
      ordered_query = T,
      significant = T,
      evcodes = T
    )
    
    if (!is.null(g)) {
      g_tbl <- g$result %>%
        filter(term_size < term_size_limit) %>%
        select(query, source, term_id, term_name, term_size, intersection, p_value) %>%
        pivot_wider(
          id_cols = c(source, term_id, term_name, term_size),
          names_from = query,
          values_from = c(intersection, p_value),
          names_vary = 'slowest'
        )
      
      g$result <- g$result %>%
        mutate(
          query = recode(
            query,
            'e' = e_title,
            's' = s_title
          )
        ) %>% 
        filter(term_size < term_size_limit)
      
      out_plot <- gostplot(g, capped = F, interactive = T)
    } else {
      out_plot <- div(style = 'height:300px;', 'No data to show.')
    }
  } else {
    out_plot <- div(style = 'height:300px;', 'No data to show.')
  }
  
  # gProfiler plot
  sketch <- htmltools::withTags(table(
    class = 'display',
    thead(
      tr(
        th(rowspan = 2, 'Source'),
        th(rowspan = 2, 'Term ID'),
        th(rowspan = 2, 'Term name'),
        th(rowspan = 2, 'Term size'),
        th(colspan = 2, style = 'background-color:#e1e8e3;border-bottom:none;text-align: center;', str_c(e_title, ' (N=', length(e_genes), ')')),
        th(colspan = 2, style = 'background-color:#e8e1e1;border-bottom:none;text-align: center;', str_c(s_title, ' (N=', length(s_genes), ')'))
      ),
      tr(
        lapply(rep(c('Intersection', 'Pval'), 2), th)
      )
    )
  ))
  
  render_sort <- "
    function(data, type, row, meta) {
      if (type === 'sort' && data === null) {
        return meta.settings.aaSorting[0][1] === 'asc' ? 999 : -999;
      }
      return data === null ? data : data.toPrecision(3);
    }
  "
  
  render_tooltip <- "
    function(row, data) {
      // 5th & 7th columns (intersections)
      for (let i = 4; i < 7; i+=2) {
        let genes = data[i] ? data[i].split('|')[1].split(',') : [];
        $('td:eq(' + i + ')', row).
          html('<span class=\"tbl-info\" title=\"' + genes.join('&#010;') + '\">' + genes.length + '</span>');
      }
    }
  "
  
  if (exists('g') && !is.null(g) && nrow(g_tbl) > 0) {
    for (query in c('e', 's')) {
      for (colname in c('intersection', 'p_value')) {
        field <- str_c(colname, query, sep = '_')
        if (!field %in% names(g_tbl)) g_tbl[[field]] <- NA
      }
    }
    
    out_tbl <- g_tbl %>%
      mutate(
        source = factor(source),
        term_id = if_else(str_detect(term_id, '^GO:'), str_c('<a href="https://www.ebi.ac.uk/QuickGO/term/', term_id, '">', term_id, '</a>'), term_id),
        intersection_e = add_count(intersection_e),
        intersection_s = add_count(intersection_s)
      ) %>%
      select(source, term_id, term_name, term_size, intersection_e, p_value_e, intersection_s, p_value_s) %>% 
      datatable(
        height = '300px',
        rownames = F,
        filter = 'top',
        container = sketch,
        class = 'compact',
        escape = F,
        options = list(
          dom = 'tip',
          scrollX = T,
          scrollY = T,
          rowCallback = JS(render_tooltip),
          columnDefs = list(
            list(targets = c(4, 6), className = 'dt-right'),
            list(targets = c(5, 7), render = JS(render_sort))
          )
        )
      ) %>%
      formatRound('term_size', 0)
  } else {
    out_tbl <- div(style = 'height:300px;', 'No data to show.')
  }
  
  list(
    tabPanel(str_c('gProfiler plot (', gp_a, ')'), out_plot),
    tabPanel(str_c('gProfiler table (', gp_a, ')'), out_tbl)
  )
}

get_tab_box <- function(typeout, cluster, genes_files, results) {
  print(str_c('get_tab_box() for ', cluster, '...'))
  
  # Biodomain modules
  imgs <- str_c(
    file.path(img_dir, typeout, round_num), 
    ifelse(is_sc, str_c(nameset, '.', cluster), nameset),
    c('TREAT-AD', 'TMT-AD', 'Mostafavi_etal', 'Milind_etal', 'Wan_etal'), 
    'Modules_Up-Down.full.logfdr.png',
    sep = '.'
  )
  names(imgs) <- c('Treat-AD', 'Tmt-AD', 'Mostafavi, et al.', 'Milind, et al.', 'Wan, et al.')
  
  biodomain_modules <- lapply(names(imgs), function(i) { tabPanel(i, img(src = imgs[[i]])) })
  
  # Biodomain correlation
  img_corr <- str_c(
    file.path(img_dir, typeout, round_num), 
    ifelse(is_sc, str_c(nameset, '.', cluster), nameset),
    'TREAT-AD.correlations.png',
    sep = '.'
  )
  
  biodomain_correlation <- list(tabPanel('Treat-AD correlation', img(src = img_corr)))
  
  # L2FC correlation
  drug_df <- genes_files[[cluster]][['drug']]
  geno_df <- genes_files[[cluster]][['geno']]
  
  pval_drug <- str_c(pval_col, 'drug', sep = '.')
  pval_geno <- str_c(pval_col, 'geno', sep = '.')
  lfc_unshrunken_drug <- str_c(lfc_unshrunken_col, 'drug', sep = '.')
  lfc_unshrunken_geno <- str_c(lfc_unshrunken_col, 'geno', sep = '.')
  
  both <- merge(drug_df, geno_df, by = gene_col, suffix = c('.drug', '.geno'))
  both <- merge(ensembl_genes, both, by.x = ensembl_var, by.y = gene_col)
  both <- both[(get(pval_drug) < 0.05 & get(pval_geno) < 0.05)]  # not adjusted p-val
  
  res <- cor.test(both[[lfc_unshrunken_drug]], both[[lfc_unshrunken_geno]])
  
  sub <- sprintf('R = %.3f; p < %.3g', res$estimate[[1]], res$p.value)
  gg_fig <- both %>%
    ggplot(aes(x = get(lfc_unshrunken_drug), y = get(lfc_unshrunken_geno))) +
    geom_point(color = NA) +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'gray') +
    geom_vline(xintercept = 0, linetype = 'dashed', col = 'gray') +
    geom_ribbon(stat = 'smooth', method = 'lm', alpha = 0.5, fill = '#fee0b6') +
    geom_smooth(method = 'lm', se = F, color = '#053061') +
    xlab(str_c(lfc_unshrunken_col, drug, sep = '.')) +
    ylab(str_c(lfc_unshrunken_col, geno, sep = '.')) +
    theme_classic()
  
  both_sig <- both[(abs(get(lfc_unshrunken_drug)) > lfc_cutoff_drug | abs(get(lfc_unshrunken_geno)) > lfc_cutoff_geno)]
  
  gg_fig <- ggplotly(gg_fig, tooltip = 'none') %>%
    add_markers(
      x = both[[lfc_unshrunken_drug]],
      y = both[[lfc_unshrunken_geno]],
      color = I('#8073ac'),
      hovertemplate = str_c(
        '<br><b>', both$external_gene_name, '</b>',
        '<br>', drug, ' L2FC: ', round(both[[lfc_unshrunken_drug]], 3), ' (p=', signif(both[[pval_drug]], 3), ')',
        '<br>', geno, ' L2FC: ', round(both[[lfc_unshrunken_geno]], 3), ' (p=', signif(both[[pval_geno]], 3), ')',
        '<extra></extra>'
      )
    )
  
  if (nrow(both_sig) > 0) {
    gg_fig <- gg_fig %>%
      add_text(
        x = both_sig[[lfc_unshrunken_drug]],
        y = both_sig[[lfc_unshrunken_geno]],
        text = both_sig$external_gene_name,
        textfont = list(size = 10),
        textposition = 'top middle'
      ) %>% add_annotations(
        x = 0,
        y = 0,
        xanchor = 'left',
        xref = 'paper',
        yref = 'paper',
        text = sub,
        showarrow = F,
        font = list(color = 'gray', size = 15)
      )
  }
  
  l2fc_corr <- list(tabPanel('L2FC correlation', gg_fig))
  
  
  # gProfiler2
  dir_enh <- results[[cluster]][category == 'enhancement' & direct == T][order(-(abs(get(lfc_geno)) + abs(get(lfc_drug))))]
  dir_sup <- results[[cluster]][category == 'suppression' & direct == T][order(-(abs(get(lfc_geno)) + abs(get(lfc_drug))))]
  
  indir_enh <- results[[cluster]][category == 'enhancement' & direct == F][order(-(abs(get(lfc_geno)) + abs(get(lfc_geno_drug))))]
  indir_sup <- results[[cluster]][category == 'suppression' & direct == F][order(-(abs(get(lfc_geno)) + abs(get(lfc_geno_drug))))]
  
  enh <- bind_rows(dir_enh, indir_enh)
  sup <- bind_rows(dir_sup, indir_sup)
  
  comp_enh <- results[[cluster]][category == 'compensatory enhancement'][order(-(abs(get(lfc_geno)) + abs(get(lfc_drug))))]
  comp_sup <- results[[cluster]][category == 'compensatory suppression'][order(-(abs(get(lfc_geno)) + abs(get(lfc_drug))))]
  
  gp_enhs <- list(enh, dir_enh, comp_enh)
  gp_sups <- list(sup, dir_sup, comp_sup)
  gp_names <- c('', 'direct', 'compensatory')
  gp_abbrevs <- c('indir + dir', 'dir', 'comp')
  
  gprofilers <- lapply(
    seq(gp_names),
    function(i) {
      get_gprofiler(gp_enhs[[i]], gp_sups[[i]], gp_names[[i]], gp_abbrevs[[i]])
    }
  ) %>% flatten()
  
  # tabBox object
  m <- do.call(tabBox, c(biodomain_modules, biodomain_correlation, l2fc_corr, gprofilers))
  m$attribs$class <- 'col-sm-12'
  
  m
}


# -----------------
# Create dashboard
# -----------------

ui <- dashboardPage(
  skin = 'black',
  header = dashboardHeader(
    title = 'DE dashboards'
  ),
  sidebar = dashboardSidebar(
    div(
      selectizeInput('cell_type', label = 'Cell type:', choices = cell_types),
    ),
    uiOutput('page_tab')
  ),
  body = dashboardBody(
    tags$head(
      includeCSS('www/assets/css/style.css')
    ),
    uiOutput('page_content')
  )
)

server <- function(input, output, session) {
  
  # Update w/ cell type selection
  clusters <- reactive({
    get_clusters(input$cell_type)
  }) %>%
    bindCache(input$cell_type)
  
  genes_files <- reactive({
    get_genes_files(input$cell_type, clusters())
  }) %>%
    bindCache(str_c(input$cell_type, '_genes'))
  
  results <- reactive({
    get_results(clusters(), genes_files())
  }) %>%
    bindCache(str_c(input$cell_type, '_res'))

  # Render sidebar tabs
  output$page_tab <- renderUI({
    cls <- names(clusters())
    subclass <- cls[[1]]
    
    pages <- lapply(
      cls,
      function(cl) {
        m <- menuItem(
          cl, 
          tabName = str_c('page_', cl), 
          icon = icon(if_else(cl == subclass, 'brain', 'dna'))
        )
        
        if (cl == subclass) {
          m <- m %>% tagAppendAttributes(class = 'active')
        }
        
        m
      }
    )
    
    do.call(sidebarMenu, pages)
  }) %>% 
    bindCache(str_c(input$cell_type, '_tab'))
  
  # Render tabs content
  output$page_content <- renderUI({
    res <- results()
    cls <- names(res)
    subclass <- cls[[1]]

    pages <- lapply(
      cls,
      function(cl) {
        tabItem(
          tabName = str_c('page_', cl),
          class = if_else(cl == subclass, 'active', ''),
          fluidRow(
            get_tab_box(subclass, cl, genes_files(), results())
          ),
          fluidRow(
            tabBox(
              width = 12,
              tabPanel(
                'Genes list',
                res[[cl]] %>% 
                  mutate(
                    ensembl_gene_id = str_c('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=', ensembl_gene_id, '">', ensembl_gene_id, '</a>'),
                    external_gene_name = str_c(external_gene_name, gene_biotype, chromosome_name, sep = '|'),
                    category = factor(category),
                    treatAD = add_count(treatAD),
                    tmtAD = add_count(tmtAD),
                    mostafavi = add_count(mostafavi),
                    milind = add_count(milind),
                    wan = add_count(wan),
                    Family = str_c(Family, 'all'),
                    chembl = str_c(add_count(chembl_protein_labels, '\\|', '~'), if_else(is_mouse_chembl, 'Mouse', 'Human'), sep = '~')
                  ) %>% 
                  dplyr::select(-gene_biotype, -chromosome_name, -is_mouse_chembl, -chembl_protein_labels) %>% 
                  datatable(
                    selection = 'none',
                    rownames = F, 
                    filter = 'top',
                    container = sketch,
                    class = 'compact',
                    escape = F,
                    callback = JS("$(\"[data-type='number'] input[type='search']\").attr('placeholder','');"),
                    options = list(
                      dom = 'tip',
                      scrollX = T,
                      scrollY = T,
                      rowCallback = JS(render_tooltip),
                      columnDefs = list(
                        list(targets = modules_col:(modules_col + 6), className = 'dt-right')
                      )
                    )
                  ) %>% 
                  formatRound(grep(str_c(lfc_col, '|log2'), names(res[[cl]]), value = T), 3) %>%
                  formatSignif(grep(str_c(p_col, '|', pval_col, '|adj|val'), names(res[[cl]]), value = T), 3) %>% 
                  renderDataTable()
              ),
              tabPanel(
                'Summary table',
                res[[cl]][order(-category, -direct), .N, by = c('category', 'direct')] %>%
                  datatable(
                    width = '400px',
                    selection = 'none',
                    rownames = F,
                    colnames = c('Category', 'Direct', 'N'),
                    class = 'compact display',
                    options = list(
                      dom = 't'
                    )
                  )
              )
            )
          )
        )
      }
    )
    
    do.call(tabItems, pages)
  }) %>% 
    bindCache(str_c(input$cell_type, '_content'))
}

shinyApp(ui, server)
