# maeplyr

`maeplyr` is an R package providing `dplyr`-like tools for manipulation of 
[`MultiAssayExperiment`](http://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html)
and [`SummarizedExperiment`](http://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
objects.

## How to install

You can install the development version from GitHub with:

```
# install.packages("devtools")
devtools::install_github("MTurner0/maeplyr")
```

`maeplyr` requires the Bioconductor packages `MultiAssayExperiment` and `SummarizedExperiment` (imported by `MultiAssayExperiment`).
If `maeplyr` fails to install, try:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MultiAssayExperiment")
```

And then install `maeplyr` from GitHub as shown above.

## Support

We welcome questions and comments about the package either through
[GitHub](https://github.com/MTurner0/maeplyr/issues) or
[email](mailto:margaret.turner1416@gmail.com).
