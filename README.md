# FDR-corrected SCCA R scripts

The R scripts `with_piece-wise_constant_Sigma.R` and `with_piece-wise_constant_Sigma_and_high-dimensional.R` perform simulations applying the FDR-corrected sparse CCA procedure, described in [Gossmann et. al. *FDR-Corrected Sparse Canonical Correlation Analysis with Applications to Imaging Genomics* (2017)](https://arxiv.org/abs/1705.04312), to Gaussian data with a piece-wise constant correlation structure (under many different parameter choices).

Example simulations can be found in the `examples/` directory. Before you can run the examples you need to perform the following steps in R:

1. Install the `devtools` R package (if you don't have it installed already):
    ```R
    install.packages("devtools")
    ```

2. Open the file `FDRcorrectedSCCA.Rproj` in [RStudio](https://www.rstudio.com/).

3. Then run
    ```R
    devtools::document()
    ```

4. Still in RStudio press CTRL-Shift-B. (If any of steps 1-4 fail, you are probably missing some other R packages, see the Error or Warning messages.)

5. Run the R scripts from the `examples/` directory (either directly from your R session, or from the command line using `Rscript with_piece-wise_constant_Sigma.R` for example).
