# How causal perspectives can inform problems in computational neuroscience

## Content

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Results and Figure Reproduction](#results-and-figure-reproduction)
- [Issues](https://github.com/ebridge2/causal_neuro/issues)
- [Citation](#citation)

## Overview

Over the past two decades, considerable strides have been made in advancing neuroscience techniques, yet the translation of these advancements into clinically relevant insights for human mental health remains a challenge. This review addresses a fundamental issue in neuroscience -- attributing causality -- and advocates for the development of robust causal frameworks. We systematically introduce the necessary definitions and concepts, emphasizing the implicit role of causal frameworks in neuroscience investigations. We illustrate how persistent challenges in neuroscience, such as batch effects and selection biases, can be conceptualized and approached using causal frameworks. Through theoretical development and real-world examples, we show how these causal perspectives highlight numerous shortcomings of existing data collection strategies and analytical approaches. We demonstrate how causal frameworks can inform both experimental design and analysis, particularly for observational studies where traditional randomization is infeasible. Using neuroimaging as a detailed case study, we explore the advantages, shortcomings, and implications for generalizability that these perspectives afford to existing and novel research paradigms. Together, we believe that this perspective offers a framework for conceptualizing, framing, and inspiring innovative approaches to problems in neuroscience.

## Repo Contents

- [R](./R): R helper functions for experiments.
- [Experiments](./experiments): analysis, inference, and plotting code for the figures in our manuscript.
- [data](./data): data from simulations used during experiments.
- [Figures](./figures): figures for our manuscript.

## System Requirements

This manuscript requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB  
CPU: 4+ cores, 3+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 4 cores@3 GHz) and internet of speed 100 Mbps.

## Installation Guide

The functions contained herein were executed using `R` version `4.3.2`. With `R` installed, you can install the dependency packages from an `R` terminal (e.g., command line, or Rstudio) using:

```
packages <- c("causalBatch", "collapse", "copula", "drtmle", "EValue", 
              "ggExtra", "ggpubr", "ggridges", "gridExtra", "jsonlite",
              "lmtest", "MASS", "MatchIt", "nnet", "parallel", "patchwork",
              "sandwich", "scales", "SuperLearner", "survey", "survival",
              "tidyverse", "transport", "twang", "WeightIt")

install.packages(packages)
```

## Results and Figure Reproduction

All of the figures for our paper can be reproduced by identifying the appropriate Figures from [Experiments](https://github.com/ebridge2/causal_neuro/tree/main/Experiments). 

## Citation

To cite conclusions from this work, you can use the following MLA citation:

Bridgeford, Eric W., Brian Caffo, Maya B. Mathur, and Russell A. Poldrack. "How Causal Perspectives Can Inform Problems in Computational Neuroscience."
