# R package: tradeSeq
## TRAjectory Differential Expression analysis for SEQuencing data

tradeSeq provides a flexible method for finding genes that are differentially expressed along one or multiple trajectories, using a variety of tests suited to answer many questions of interest.

## Installation

To get the development version, run 
```r
devtools::install_github("statOmics/tradeR/")
```

## Issues and bug reports

Please use https://github.com/statOmics/tradeSeq/issues to submit issues, bug reports, and comments.

## Usage 

Start with the vignette either [online](https://statOmics.github.io/tradeSeq/) or by running
```{r}
utils::vignette(topic = "tradeSeq", package = "tradeSeq")
```

## Cheatsheet

You can also refer to this cheatsheet to undersand a common workflow

![](vignettes/cheatsheet_highRes.jpeg)

## Contributing and requesting

A number of tests have been implemented in tradeSeq, but researchers may be interested in other hypotheses that current implementations may not be able to address. We therefore welcome contributions on GitHub on novel tests based on the tradeSeq model.
Similar, you may also request novel tests to be implemented in tradeSeq by the developers, preferably by adding an issue on the GitHub repository. If we feel that the suggested test is widely applicable, we will implement it in tradeSeq.
