# dnam-model-power 

Calculating power to detect DNAm-trait associations in EWAS when the associations arise because of three different reasons: confounding, forward-cause (DNAm changes cause trait variation), reverse-cause (trait variation causes DNAm changes). 

## Reports 

* [PDF report](supplementary-note-EWAS-power.pdf)
* [HTML report](dnam_model_power.nb.html)

## Instructions on making reports

RMarkdown script R package dependencies:

``` r
# In R
install.packages(c("pwr", "dplyr", "ggplot2", "bookdown", "rmarkdown", "tinytex"))
## STEP BELOW IS ONLY NEEDED TO PRODUCE THE PDF DOCUMENT
tinytex::install_tinytex()
```

You also require [`preamble.tex`](preamble.tex) to create the PDF.

To make the original html doc Gib wrote:

``` bash
Rscript -e "rmarkdown::render('dnam_model_power.rmd', output_format = 'all')"
```

To make the paper-ready supplementary note version:

``` bash
Rscript -e "rmarkdown::render('dnam_model_power_pdfout.rmd', output_format = 'all', output_file = 'supplementary-note-EWAS-power.pdf')"
```