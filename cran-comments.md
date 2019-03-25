## Test environments

* local Ubuntu 16.04 LTS install, R 3.5.3
* Ubuntu 16.04 LTS, R-devel (via rhub::check)
* devtools::check_win_release

## R CMD check results

### Local Ubuntu 16.04 LTS Install (R 3.5.3)

Command `R CMD check fastGeneMI_1.0.tar.gz --as-cran`  output:

There were no ERRORs. 

There was 1 WARNING:

* checking compilation flags used ... WARNING
Compilation used the following non-portable flag(s):
  ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’

This is due to the local compiler and does not affect CRAN.

There was 1 NOTE:

Package has a FOSS license but eventually depends on the following
package which may restrict use:
  minet

The [minet package](http://www.bioconductor.org/packages/release/bioc/html/minet.html) has
"file LICENSE" in the License field on Bioconductor. This may be a typo in the minet
package.

### Ubuntu 16.04 LTS (R-devel, via rhub::check)

Estimator tests fail.

### devtools::check_win_release

No ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Jonathan Ish-Horowicz <jonathan.ish-horowicz17@imperial.ac.uk>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  ARACNE (10:59)
  CLR (10:54)
  MRNET (10:70)
  Rcpp (11:43)
  RcppArmadillo (11:52)

Package has a FOSS license but eventually depends on the following
package which may restrict use:
  minet


Those words are spelt correctly.

The [minet package](http://www.bioconductor.org/packages/release/bioc/html/minet.html) has
"file LICENSE" in the License field on Bioconductor. This may be a typo in the minet
package.

## Downstream dependencies

There are currently no downstream dependencies for this package.