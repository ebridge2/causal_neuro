FROM bioconductor/bioconductor_docker:devel

RUN Rscript -e "install.packages(c('tidyverse', 'parallelDist', 'MatchIt', 'energy', 'multcomp', 'survey', 'dplyr'))"
RUN Rscript -e "install.packages(c('WeightIt', 'SuperLearner', 'drtmle', 'causalBatch', 'copula'))"
