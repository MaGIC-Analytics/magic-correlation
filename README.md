# MaGIC Gene Correlation & Co-expression Tool

![GitHub last commit](https://img.shields.io/github/last-commit/MaGIC-Analytics/magic-correlation)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![made with Shiny](https://img.shields.io/badge/R-Shiny-blue)](https://shiny.rstudio.com/)

Compute and visualize gene-gene correlations from an expression matrix. Supports a single-gene-vs-all-genes mode with analytical p-values and BH-adjustment, plus an N-gene correlation matrix mode with ComplexHeatmap and corrplot visualizations.

## Running the App
This Shiny App has been built in to a docker container for easy deployment. You can build the image yourself (and thereby customize any ports you need) after downloading it:
```
docker build -t correlation .
docker run -d --rm -p 8080:8080 correlation
#Or for testing docker run -t -i --rm -p 8080:8080 correlation
```
And it should be hosted at localhost:8080
