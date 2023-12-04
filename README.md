# Introduction to latent means models in R

Repository for materials supporting our manuscript discussing practical
implementations of latent mean centering in multilevel models of
longitudinal psychological data.

- Repo on GitHub: <https://github.com/mvuorre/latent-mean-centering-ms>
- Archived repo: tbd
- Online: tbd

## Reproduce

Everything is coded in [R](https://cran.r-project.org/) and bunched
together as a [Quarto](https://quarto.org/) website. To reproduce, clone
the repo and render the Quarto documents:

``` bash
# Download materials
git clone https://github.com/mvuorre/latent-mean-centering-ms.git
cd latent-mean-centering-ms

# Set up R environment
Rscript -e 'renv::restore(prompt = FALSE)'

# Run code and render outputs
quarto render
```

## Contribute

Contributions preferably via pull requests.
