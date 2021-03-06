# Hickory

### Bayesian estimates of Wright's F statistics

This package estimates the within population inbreeding coefficient,
<img src="https://render.githubusercontent.com/render/math?math=f">,
and Wright's <img
src="https://render.githubusercontent.com/render/math?math=F_{ST}">,
<img
src="https://render.githubusercontent.com/render/math?math=\theta">,
using a Bayesian approach. Right now, it is limited to two alleles per
locus for co-dominant markers. It's easiest to use if your data are in
a CSV file, but it wouldn't be too difficult to convert data from
other formats to the one that required by `analyze_codominant()` and
`analyze_dominant()`. If you look at the _Roadmap_ below, you'll see
that I plan to add an interface to `adegenet`, which should allow you
to work with data in most of the widely used formats.

### A note on requirements

You will need to have a C/C++ compiler installed in order to use this
package (for now).<sup>1</sup> The
[RStan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) page has helpful information on getting your system
configured. I'll help if I can, but I've been using Macs for nearly 10
years. I can probably help with Mac or Linux problems, but I may not
be much help with Windows. If you've managed to install other R
packages from source, you'll probably be just fine.

I also specified minimum versions for all of the libraries that `Hickory`
depends on. It's possible I could get away with earlier versions, but
I don't have a way to test that. That means you might need to upgrade
some of your libraries before you can install this one.

### Installation

If you don't already have the `devtools` package installed, first
you'll need to install it.

```
install.packages("devtools")
```

Once that's installed, then you can install `Hickory` like this:

```
install.packages(c("bayesplot", "rstan", "tidyverse"))
devtools::install_github("kholsinger/Hickory")
```

If you also want to install vignettes illustrating how to use
`Hickory` in more detail, you'll want to change that second line to

```
devtools::install_github("kholsinger/Hickory", build_vignettes = TRUE)
```

The installation will take longer, but you'll be able to run
`vignette("Hickory")` to get an overview, and
`browseVignettes("Hickory")` to see all of the vignettes that are
available. If you'd prefer, you can simply look at the vignettes in
the `doc` directory here: 

[https://github.com/kholsinger/Hickory/tree/master/doc](https://github.com/kholsinger/Hickory/tree/master/doc)

You'll probably be interested only in the HTML files in this
directory. To view them, you'll need to hit the "Raw" button on the
upper right of where the HTML code is displayed, save the file on you
computer, and open it from your computer.

The Rmd file is the [rmarkdown](https://rmarkdown.rstudio.com/) file
that produced the HTML. The R file is R code that's extracted from
the Rmd file.

### A note on dominant markers

Foll et al.<sup>2</sup> identified biases associated with the method
originally implemented in `Hickory` for dominant markers. Several
users also reported that estimates of f seemed unreasonably high when
they used a large number of markers. 

The C++ version of `Hickory` used Metropolis-Hastings sampling to
approximate the posterior. Sampling was slow and estimates of <img
src="https://render.githubusercontent.com/render/math?math=\theta">
showed high autocorrelation. The Hamiltonian Monte Carlo algorithm
used in Stan doesn't suffer from either of those limitations, but I
haven't checked the results from this sampler against the biases that
Foll et al. reported. Use caution interpreting estimates of <img
src="https://render.githubusercontent.com/render/math?math=f"> until
I've checked that out. Fortunately, estimates of <img
src="https://render.githubusercontent.com/render/math?math=\theta">
aren't too sensitive to <img
src="https://render.githubusercontent.com/render/math?math=f">, so
those estimates are likely to be reliable.

If you plan to use `Hickory` to analyze dominant markers, I encourage
you to read the vignette, `vignette("dominant")`.

### A note on the Github repository

If you visit the repository you'll see that the "main" branch is still
called "master". I plan to rename it after Github releases tools
making it easier. I am not good at `git`, and I need all of the help I
can get.

### Roadmap

There are a lot of improvements that still need to be made. In
addition to the items listed below, I'll be working on improved
documentation. Please don't hesitate to email me if you have a
question, a comment, or suggestions for future improvements.

1. Implement <img
   src="https://render.githubusercontent.com/render/math?math=f=0">
   and 
   <img
   src="https://render.githubusercontent.com/render/math?math=\theta=0">
   models with model comparison using 
   `loo` to provide ways to evaluate whether there is evidence for
   inbreeding and whether there is evidence for allele frequency
   differences among populations. 
   
   ***DONE***
   
2. Implement population- and locus-specific effects on <img
   src="https://render.githubusercontent.com/render/math?math=\theta">,
   including identification of potential outlier loci and populations.
   
   ***DONE***

3. Implement interface with `adegenet`.

4. Build in some internal tests to ensure accuracy and consistency.

5. Implement posterior predictive checks.

6. Investigate possible biases in dominant marker estimates.
   
7. Implement multiallele version of `analyze_codominant()`.

8. Implement locus and population selection

<sup>1</sup>Once I feel comfortable enough with this package to release it
    to CRAN, I believe that the build system there will produce the
    binaries for different platforms. I don't have the ability to do
    that myself.
    
<sup>2</sup>Foll M, Beaumont MA, Gaggiotti O. 2008. An Approximate
    Bayesian Computation Approach to Overcome Biases That Arise When
    Using Amplified Fragment Length Polymorphism Markers to Study
    Population Structure. <em>Genetics</em> 179:927-939.
