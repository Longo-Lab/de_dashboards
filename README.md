# DE dashboards

## Options

Use `rmarkdown::render()` to knit the `dashboard.Rmd` template with the following options:

- **`base_dir`**: Path to project directory
  - Must contain subfolder(s) named by `typeout` for single-cell data or `nameset` for bulk data
  - Each subfolder must contain the respective DE results files named as: `R#.nameset.(typeout).RNA.analyses.csv.gz`
- `templates_dir`: Path to templates directory (default: [`templates`](https://github.com/Longo-Lab/de_dashboards/tree/main/templates))
- **`ensembl_genes`**: Data table containing Ensembl genes + modules info
  - Can be obtained by calling `get_genes_info()` from [`genes_info.R`](https://github.com/Longo-Lab/de_dashboards/blob/main/scripts/genes_info.R)
- **`nameset`**: Name of project (e.g., `PS19_C31`)
- **`geno`**: Name of genotype (e.g., `PS19`)
- **`drug`**: Name of drug (e.g., `C31`)
- `round_num`: Round of source file (default: `R1`)
- `typeout`: Cell type if single-cell data
  - If unset, `nameset` is used
- `analyses`: Name of comparisons (default: `c('Tg-VvsWt-V', 'Tg-DvsWt-V', 'Tg-DvsTg-V')`)
  - Provide as a vector of comparisons in the following order: geno, geno + drug (stim), drug (stim)
  - Optionally, include a 4th comparison for drug (stim) effect in wildtype
  - Specifies how the DE results files are named (see above)
- `de_cols`: Column names for gene ID, log fold change, p-value, followed by any additional columns of interest (default: `c('GENE', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2')`)
  - Log fold change will be used for ranking genes in gProfiler input and determining gene category
  - P-value will be used for defining significance when determining gene category
- `de_names`: Column labels for `de_cols` not including gene ID (default: `c('Log2FC', 'Pval (adj)', 'Pct 1', 'Pct 2')`)
- `pval_col`: Column name for nominal p-value (default: `p_val`)
  - Used in the `L2FC correlation` tab
- `lfc_unshrunken_col`: Column name for unshrunken log fold change (default: `avg_log2FC`)
  - Used in the `L2FC correlation` tab
- `fdr_cutoff`: Significance value threshold (default: `0.05`)
  - Used for p-value specified in `de_cols` to determine gene category
- `lfc_cutoff_geno`: Absolute genotype log fold change cutoff (default: `0.55`)
  - Used for displaying gene label in the `L2FC correlation` tab
- `lfc_cutoff_drug`: Absolute drug log fold change cutoff (default: `0.5`)
  - Used for displaying gene label in the `L2FC correlation` tab

**Note**: Required options are in **bold**.

## Notes

Load R and run script to generate dashboards:

```
$ ml R/4.0
$ Rscript scripts/PS19_C31.R
```

Ensure [`tidyr`](https://tidyr.tidyverse.org/news/index.html#tidyr-120) is updated to at least v.1.2.0 to use the `names_vary` argument for `pivot_wider()`:

```
$ R
> install.packages('tidyr')
```

If running into [`ReferenceError: FlexDashboard is not defined`](https://github.com/rstudio/flexdashboard/issues/116), make sure `pandoc` is up-to-date. This can be done by specifying the [`RSTUDIO_PANDOC`](https://stackoverflow.com/a/29710643/6373540) environment variable in `~/.bash_profile`:

```sh
export RSTUDIO_PANDOC=/usr/lib/rstudio-server/bin/quarto/bin
```

To encrypt the HTML files using [`staticrypt`](https://www.npmjs.com/package/staticrypt), install the CLI locally with [`npm`](https://stackoverflow.com/a/14469516/6373540):

```
$ mkdir -p .npm/node_modules
$ npm install --prefix .npm/ staticrypt
```

Make sure to append to the `PATH` environment variable in `~/.bash_profile`:

```sh
PATH=$PATH:$HOME/.npm/node_modules/.bin
export PATH
```

Encrypt all HTML files in the output directory, making sure to purge modules if there are conflicts:

```
$ module purge
$ find outputs/PS19_C31/ -type f -name "*.html" -exec staticrypt {} PASSPHRASE -o {} -r 1 \;
```

Alternatively, run `gen_dashboards.sh` and provide the name (`-n`) and password (`-p`) arguments:

```
$ sbatch gen_dashboards.sh -n PS19_C31 -p PASSPHRASE
```

## Save dashboard file

| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | biodomain_correlation.R | biodomain_enrichment.R | save_dashboard_files.R | save_stim_dashboard_files.R |
|---|---|---|---|---|
| `-n`<br>`--nameset` | Name of project | Name of project | Name of project | Name of project |
| `-g`<br>`--genotype` | Name of genotype | Name of genotype | Name of genotype | Name of genotype |
| `-d`<br>`--drug` | Name of drug | Name of drug | Name of drug | Name of drug |
| `-r`<br>`--round-number` | Round of source file<br><br>**default**: `R1` | Round of source file<br><br>**default**: `R1` | Round of source file<br><br>**default**: `R1` | Round of source file<br><br>**default**: `R1` |
| `-t`<br>`--type` | Cell type if single-cell data, otherwise `bulk` for bulk data or `stim` for stimulation analysis<br><br>**default**: `bulk` | Cell type if single-cell data, otherwise `bulk` for bulk data or `stim` for stimulation analysis<br><br>**default**: `bulk` | Cell type if single-cell data, otherwise `bulk` for bulk data<br><br>**default**: `bulk` | -- |
| `-a`<br>`--analyses` | Name of comparisons (_ORDER BY geno, drug, geno + drug OR stim in wtv, tgv, tgd_)<br><br>If `-t` is NOT `stim`, **default**: `Tg-VvsWt-V,Tg-DvsTg-V,Tg-DvsWt-V`<br><br>If `-t` is `stim`, **default**: `Wt-StvsWt-Un,Tg-StvsTg-Un,Tg-D-StvsTg-D-Un` | Name of comparisons (_ORDER BY geno, drug, geno + drug OR stim in wtv, tgv, tgd_)<br><br>If `-t` is NOT `stim`, **default**: `Tg-VvsWt-V,Tg-DvsTg-V,Tg-DvsWt-V`<br><br>If `-t` is `stim`, **default**: `Wt-StvsWt-Un,Tg-StvsTg-Un,Tg-D-StvsTg-D-Un` | Name of comparisons (ORDER BY geno, drug, geno + drug)<br><br>**default**: `Tg-VvsWt-V,Tg-DvsTg-V,Tg-DvsWt-V` | Name of comparisons (ORDER BY stim in wtv, tgv, tgd)<br><br>**default**: `Wt-StvsWt-Un,Tg-StvsTg-Un,Tg-D-StvsTg-D-Un` |
| `-s`<br>`--shrunken` | Flag indicating shrunken LFC is being used<br><br>**default**: `FALSE`<br><br>_**Note**: This is more for consistency w/ the other scripts and to use default `-c` value of `shrunkenL2FC`. If providing custom shrunken LFC column name, including `-s` or not technically does not make a difference._ | Flag indicating shrunken LFC is being used<br><br>**default**: `FALSE` | Flag indicating shrunken LFC is being used<br><br>**default**: `FALSE` | Flag indicating shrunken LFC is being used<br><br>**default**: `FALSE` |
| `-u`<br>`--unadjusted-p` | -- | Flag indicating unadjusted (nominal) p-value is being used for filtering instead of adjusted p-value (in the case of unshrunken LFC)<br><br>**default**: `FALSE`<br><br>_**Note**: This is more for consistency w/ **save_dashboard_files.R** and to use default `-c` value for p-value filter. If providing custom p-value column name, including `-u` or not technically does not make a difference._ | Flag indicating unadjusted (nominal) p-value is being used for filtering instead of adjusted p-value (in the case of unshrunken LFC)<br><br>**default**: `FALSE` | -- |
| `-i`<br>`--gene-col` | Column name for gene ID (bulk) or symbol (single-cell)<br><br>**default**: `GENE` | Column name for gene ID (bulk) or symbol (single-cell)<br><br>**default**: `GENE` | Column name for gene ID (bulk) or symbol (single-cell)<br><br>**default**: `GENE` | Column name for gene ID<br><br>**default**: `GENE` |
| `-c`<br>`--de-columns` | Column name for log fold change (can provide shrunken or unshrunken)<br><br>If `-s`, **default**: `shrunkenL2FC`<br><br>If `-t` is `bulk`/`stim`, **default**: `log2FoldChange`<br><br>If `-t` is single-cell, **default**: `avg_log2FC` | If shrunken (and not stim), column name for shrunken log fold change OR if unshrunken (or stim), column names for log fold change and p-value column for filtering<br><br>If `-s` and `-t` is NOT stim, **default**: `shrunkenL2FC`<br><br>If `-s` and `-t` is stim, **default**: `shrunkenL2FC,padj`<br><br>If `-t` is `bulk`/`stim`, **default**: `log2FoldChange,padj` OR `log2FoldChange,pvalue` if `-u`<br><br>If `-t` is single-cell, **default**: `avg_log2FC,p_val_adj` OR `avg_log2FC,p_val` if `-u`<br><br>_**Note**: Must use padj < 0.05 threshold w/ either shrunken or unshrunken LFC for stim dashboard._ | If shrunken, column names in the order of `sLFC,sval,LFC,pval,padj` OR if unshrunken, column names in the order of `LFC,p1,p2,(pct1,pct2)`, where `p1` is the same p-value column used for filtering in **biodomain_enrichment.R**, `p2` is the other (out of p-value and p-adj columns), and `pct1`/`pct2` only applies for single-cell data<br><br>If `-s`, **default**: `shrunkenL2FC,svalue,log2FoldChange,pvalue,padj`<br><br>If `-t` is `bulk`/`stim`, **default**: `log2FoldChange,padj,pvalue` OR `log2FoldChange,pvalue,padj` if `-u`<br><br>If `-t` is single-cell, **default**: `avg_log2FC,p_val_adj,p_val,pct.1,pct.2` OR `avg_log2FC,p_val,p_val_adj,pct.1,pct.2` if `-u` | If shrunken, column names in the order of `sLFC,padj,sval,LFC,pval` OR if unshrunken, column names in the order of `LFC,adj,pval`<br><br>If `-s`, **default**: `shrunkenL2FC,padj,svalue,log2FoldChange,pvalue`<br><br>If NOT `-s`, **default**: `log2FoldChange,padj,pvalue` |
| `-l`<br>`--lfc-threshold` | -- | LFC threshold used to filter genes<br><br>**default**: `log2(1.1)` | LFC threshold used to filter genes<br><br>**default**: `log2(1.1)` | LFC threshold used to filter genes<br><br>**default**: `log2(1.1)` |
| `-p`<br>`--p-threshold` | -- | P-value threshold used to filter genes<br><br>**default**: `0.05`<br><br>_**Note**: Do not change default for stim dashboard._ | P-value threshold used to filter genes<br><br>**default**: `0.05` | -- |

## Shiny version

Dashboard files can be saved by running [`save_dashboard_files.R`](https://github.com/Longo-Lab/scripts/blob/main/save_dashboard_files.R) in the `base_dir`. It will read in a file named `config.R` located inside that directory. Then, if there is a file named `config.R` inside the project subfolder, it will read that in after, which can override any settings in the first file.

Specifically, the following **bolded** options must be defined, whereas other options have default values:

- **`nameset`** (needed only for single-cell data)
- **`geno`**
- **`drug`**
- `round_num` (default: `R1`)
- `analyses` (default: `c('Tg-VvsWt-V', 'Tg-DvsWt-V', 'Tg-DvsTg-V')`)
- `de_cols` (default: `c('GENE', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2')`)
- `de_names` (default: `c('Log2FC', 'Pval (adj)', 'Pct 1', 'Pct 2')`)
- `pval_col` (default: `p_val`)
- `lfc_unshrunken_col` (default: `avg_log2FC`)
- `fdr_cutoff` (default: `0.05`)
- `lfc_cutoff_geno` (default: `0.55`)
- `lfc_cutoff_drug` (default: `0.5`)
- **`is_sc`** (`T` or `F` as to whether or not it is single-cell data)

**Note**: The subfolder name (i.e., `typeout` for single-cell data or `nameset` for bulk data) is passed as a command line argument to `save_dashboard_files.R`. `ensembl_ver` and gProfiler `term_size_limit` are also set at the top of the script. The Ensembl genes data is saved inside `reference/biodomains/` by the `genes_info.R` script located there.

In the `base_dir`, add a copy of [`app.R`](https://github.com/Longo-Lab/de_dashboards/blob/main/app.R) and [`www/`](https://github.com/Longo-Lab/de_dashboards/tree/main/www). At the top of `app.R`, modify the dropdown selection choices accordingly and define the following:

- **`nameset`** (needed only for single-cell data)
- **`round_num`**
- **`is_sc`**

See [`deploy.R`](https://github.com/Longo-Lab/de_dashboards/blob/main/deploy.R) for template of deploying to [shinyapps.io](https://www.shinyapps.io/). Make sure to add this file to `.gitignore` if using as it contains credential information.

## Use

This repository is for internal use of the Longo Lab only. All work herein is under the exclusive copyright of the Longo Lab. Open source versions of the code may be made available in dedicated repositories associated with individual publications and projects at the discretion of the Longo Lab.
