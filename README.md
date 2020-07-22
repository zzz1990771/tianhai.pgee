# tianhai.pgee

This package was modified upon the the PGEE package from https://github.com/cran/PGEE. This modification that aims at accelerating the estimation speed and optimize the inplementation further is designed to accomedate the computing requirement of a non-linear model estiamtion, which use linear PGEE as one of the core components. 

This package is normally faster than the original PGEE by 50% to 500% (depends on different use cases).

# Installtion of this package

If you have the devtools package already, you can install the package from github directly:

```R
devtools::install_github("kbroman/broman")
```

Otherwise, you need to install the devtools package first:

```R
install.packages("devtools")
```
