# gonovax

<!-- badges: start -->
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![R build status](https://github.com/mrc-ide/gonovax/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/gonovax/actions)
[![CodeFactor](https://www.codefactor.io/repository/github/mrc-ide/gonovax/badge?s=1d60596994e72a75be157e74ec8e23948e90fc31)](https://www.codefactor.io/repository/github/mrc-ide/gonovax)
[![codecov](https://codecov.io/gh/mrc-ide/gonovax/branch/master/graph/badge.svg?token=9u8S3v45AX)](https://codecov.io/gh/mrc-ide/gonovax)
<!-- badges: end -->

This package implements a compartmental model of gonorrhoea infection with vaccination.

<img src="man/figures/vax_diagram.png" align="center" style = "border: none; float: center;" width = "800px">


## Installation

You will need a compiler to install dependencies for the package, and to build
the models. Use `pkgbuild::check_build_tools()` to see if your system is usable.

You will need the packages `odin` and `mcstate`, which can be installed using:

```r
remotes::install_github("mrc-ide/odin", upgrade = FALSE)
remotes::install_github("mrc-ide/mcstate", upgrade = FALSE)
```


The package can then be installed directly from GitHub with:

```r
remotes::install_github("mrc-ide/gonovax", upgrade = FALSE)
```

## License

MIT © Imperial College of Science, Technology and Medicine

