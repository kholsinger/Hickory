# Hickory

### Bayesian estimates of Wright's F statistics

This package estimates the within population inbreeding coefficient,
$f$, and Wright's $F_{ST}$, $\theta$, using a Bayesian approach. Right
now, it is limited to two alleles per locus for co-dominant
markers. It's easiest to use if your data are in a CSV file, but it
wouldn't be too difficult to convert data from other formats to the
one that required by `analyze_codominant()` and
`analyze_dominant()`. If you look at the _Roadmap_ below, you'll see
that I plan to add an interface to `adegenet`, which should allow you
to work with data in most of the widely used formats.

### A note on requirements

You will need to have a C/C++ compiler installed in order to use this
package (for now).[^1] The (RStan Getting
Started)[https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started]
page has helpful information on getting your system configured. I'll
help if I can, but I've been using Macs for nearly 10 years. I can
probably help with Mac or Linux problems, but I may not be much help
with Windows. If you've managed to install other R packages from
source, you'll probably be just fine. 

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
install.packages("bayesplot", "rstan", "tidyverse")
devtools::install_github("kholsinger/Hickory")
```

### Roadmap

There are a lot of improvements that I need to make before this
package is ready. Here's what I have in mind, in rough sequence.

1. Implement $f=0$ and $\theta=0$ models with model comparison using
   `loo` to provide ways to evaluate whether there is evidence for
   inbreeding and whether there is evidence for allele frequency
   differences among populations.
2. Implement population- and locus-specific effects on $\theta$,
   including identification of potential outlier loci.
3. Implement interface with `adegenet`.
4. Implement multiallele version of `analyze_codominant()`.



[^1]: Once I feel comfortable enough with this package to release it
    to CRAN, I believe that the build system there will produce the
    binaries for different platforms. I don't have the ability to do
    that myself.
