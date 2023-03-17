# ABXbug

This repository collects the code for the analyses presented in the paper: [Unravelling the collateral damage of antibiotics on gut bacteria](https://doi.org/10.1038/s41586-021-03986-2).

The accompanying data needed to reproduce the results are available through Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3527540.svg)](https://doi.org/10.5281/zenodo.3527540)
 
## Reproduction

In order to reproduce the results presented in the paper, you can use the provided scripts. All analyses were performed using R version `4.1.0` and the needed packages can be installed with the provided `renv.lock` file, since this project uses the project environment management package [`renv`](https://rstudio.github.io/renv/index.html). You can check out the documentation for more information: [renv tutorial](https://rstudio.github.io/renv/articles/renv.html).

After installing all required packages, you can automatically download the data by executing the `setup.R` script in the `utils` folder.

# Tidy.py

+ It is used to choose the sample sequence form a list and diveided the sequences into different .txt. 
+ Use the tool.bat to change the name form `.txt` to `.fasta`.





## Contact

For any questions or comments, please contact the corresponding author or `l.maier [at] uni-tuebingen.de`.
