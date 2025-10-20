# How causal perspectives can inform problems in neuroscience data analysis

## Content

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Results and Figure Reproduction](#results-and-figure-reproduction)
- [Issues](https://github.com/ebridge2/causal_neuro/issues)
- [Citation](#citation)

## Overview

Over the past two decades, considerable strides have been made in advancing neuroscientific techniques, yet challenges remain in attributing causality to observed associations. This review addresses a fundamental issue in observational neuroscience studies and advocates for incorporating causal inference frameworks into standard practice. We systematically introduce necessary definitions and concepts, emphasizing how causal assumptions underlie statistical analyses even when not explicitly stated. Through a running example on sleep quality and white matter integrity, we illustrate how persistent challenges, including confounding and selection biases, can be conceptualized and addressed using causal frameworks. We demonstrate practical approaches for making assumption violations transparent through hands-on examples: supplementary case studies using multi-site harmonization and head motion exclusion procedures provide step-by-step diagnostic techniques for checking covariate overlap and identifying selection bias through exclusion pattern analysis. We explore how these causal perspectives can inform both experimental design and analytical choices, particularly for observational studies where traditional randomization is infeasible. Together, we believe this framework offers concrete tools for strengthening causal interpretations and inspiring more robust approaches to problems in neuroscience.

## Repo Contents

- [Exploratory](./Exploratory): Rmarkdown documents with detailed causal exploratory analyses.
- [Figures](./Figures): figures for our manuscript.

## System Requirements

This manuscript requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB  
CPU: 4+ cores, 3+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 4 cores@3 GHz) and internet of speed 100 Mbps.

## Installation Guide

The functions contained herein were executed using `R` version `4.3.2`. With `R` installed, you can install the dependency packages from an `R` terminal (e.g., command line, or Rstudio) using:

```r
packages <- c("collapse", "dagitty", "ggdag", "ggExtra", "ggridges", 
              "gridExtra", "jsonlite", "knitr", "patchwork", "tidyverse")

install.packages(packages)
```

## Results and Figure Reproduction

All of the figures for the Appendix of our paper can be reproduced by identifying the appropriate section from [Exploratory](https://github.com/ebridge2/causal_neuro/tree/main/Exploratory). For real-data analyses using the ABCD study, you will need the appropriate files indicated for the corresponding analyses placed in the folder `<path/to/this/repository>/data/abcd`.

## Citation

To cite conclusions from this work, you can use the following MLA citation:

Bridgeford, Eric W., Brian Caffo, Maya B. Mathur, and Russell A. Poldrack. "How Causal Perspectives Can Inform Problems in Neuroscience Data Analysis." arXiv (2025).
