# DE dashboards

## Notes

Load R and run script to generate dashboards:

```
$ ml R/4.0
$ Rscript ps19_c31.R
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

Encrypt all HTML files in the output directory, making sure to purge modules if there are conflicts with R:

```
$ module purge
$ find ./outputs/ps19_c31/ -type f -name "*.html" -exec staticrypt {} PASSPHRASE -o {} -r 1 \;
```

Alternatively, run `gen_dashboards.sh` and provide the name (`-u`) and password (`-p`) arguments:

```
$ sbatch gen_dashboards.sh -n ps19_c31 -p PASSPHRASE
```

## Use

This repository is for internal use of the Longo Lab only. All work herein is under the exclusive copyright of the Longo Lab. Open source versions of the code may be made available in dedicated repositories associated with individual publications and projects at the discretion of the Longo Lab.
