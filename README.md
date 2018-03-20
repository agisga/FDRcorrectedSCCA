# FDR-corrected SCCA R scripts

The R scripts in the `examples` directory perform the simulation studies described in [Gossmann et. al. *FDR-Corrected Sparse Canonical Correlation Analysis with Applications to Imaging Genomics* (2017)](https://arxiv.org/abs/1705.04312). Code to generate visualizations from the simulation results, as presented in the paper, is provided in this repository as well.

The R scripts have been written with the intention of running them on the [high performance computing cluster Cypress at Tulane University](https://crsc.tulane.edu/), which uses the [Slurm resource management system](https://slurm.schedmd.com/documentation.html).

## R package

Many functions (see the `R` directory), which are used to implement the simulations in the `examples` directory, are provided in form of an R package. The R package `FDRcorrectedSCCA` can be installed using `devtools` in R (see below). However, this should not be necessary in order to run the scripts from the `examples` directory, because within those scripts all functions from the `FDRcorrectedSCCA` package are loaded with `devtools::load_all()`. However, if you have trouble running the example scripts (*ahem* Windows *ahem*...), try installing the required functions as an R package following the instructions below.

If you wish to install the R package `FDRcorrectedSCCA`, then perform the following steps in R:

1. Install the `devtools` R package (if you don't have it installed already):
    ```R
    install.packages("devtools")
    ```

2. Open the file `FDRcorrectedSCCA.Rproj` in [RStudio](https://www.rstudio.com/) (you can create the Rproj file from RStudio via `File > New Project > Existing Directory`).

3. Then run
    ```R
    devtools::document()
    ```

4. Still in RStudio press CTRL-Shift-B. (If any of steps 1-4 fail, you are probably missing some other R packages, see the Error or Warning messages. Install the missing packages and try again from step 2.)
