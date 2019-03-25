## Test environments

* local ubuntu 16.04 install, R 3.4.4
* local Windows 10 install, R 3.4.4
* devtools::check_win_release

## R CMD check results

### Local ubuntu 16.04 Install (R 3.4.4)

Command R CMD check --as-cran output:

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

Package has a FOSS license but eventually depends on the following
package which may restrict use:
  minet

The minet package (http://www.bioconductor.org/packages/release/bioc/html/minet.html) has
"file LICENSE" in the License field on Bioconductor. This may be a typo in the minet
package.

### Local Windows 10 install (R 3.4.4)

devtools::check() output:

0 errors √ | 0 warnings √ | 0 notes √

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

The minet package (http://www.bioconductor.org/packages/release/bioc/html/minet.html) has
"file LICENSE" in the License field on Bioconductor. This may be a typo in the minet
package.

## Downstream dependencies

There are currently no downstream dependencies for this package.



