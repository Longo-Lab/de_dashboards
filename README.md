# DE dashboards

## Options

Use `rmarkdown::render()` to knit the `dashboard.Rmd` template with the following options:

- **`base_dir`**: Path to project directory
  - Must contain subfolder(s) named by `typeout` for single-cell data or `nameset` for bulk data
  - Each subfolder must contain the respective DE results files named as: `filedate.nameset.(typeout).RNA.analyses.csv.gz`
- `templates_dir`: Path to templates directory (default: [`templates`](https://github.com/Longo-Lab/de_dashboards/tree/main/templates))
- **`ensembl_genes`**: Data table containing Ensembl genes + modules info
  - Can be obtained by calling `get_genes_info()` from [`genes_info.R`](https://github.com/Longo-Lab/de_dashboards/blob/main/scripts/genes_info.R)
- **`nameset`**: Name of project (e.g., `PS19_C31`)
- **`geno`**: Name of genotype (e.g., `PS19`)
- **`drug`**: Name of drug (e.g., `C31`)
- **`filedate`**: Date of source file (e.g., `20220711`)
- `typeout`: Cell type if single-cell data
  - If unset, `nameset` is used
- `analyses`: Name of comparisons (default: `c('Tg-VvsWt-V', 'Tg-DvsWt-V', 'Tg-DvsTg-V')`)
  - Provide as a vector of comparisons in the following order: geno, geno + drug (stim), drug (stim)
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
$ Rscript scripts/ps19_c31.R
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
$ find outputs/ps19_c31/ -type f -name "*.html" -exec staticrypt {} PASSPHRASE -o {} -r 1 \;
```

Alternatively, run `gen_dashboards.sh` and provide the name (`-n`) and password (`-p`) arguments:

```
$ sbatch gen_dashboards.sh -n ps19_c31 -p PASSPHRASE
```

## Use

This repository is for internal use of the Longo Lab only. All work herein is under the exclusive copyright of the Longo Lab. Open source versions of the code may be made available in dedicated repositories associated with individual publications and projects at the discretion of the Longo Lab.
