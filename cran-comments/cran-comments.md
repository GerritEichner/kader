---
title: "cran-comments.md"
author: "Gerrit Eichner"
date: "24 August 2017"
output: html_document
---

## Test environments
* local OS X install, R 3.1.2 ???????????
* ubuntu 12.04 (on travis-ci), R 3.1.2 ?????????????
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:

* Possibly mis-spelled words in DESCRIPTION:  
  Eichner (12:16, 13:32)  
  Srihera (11:73)  
  Stute (11:63, 12:28, 13:44)  
  nonparametric (10:67)  
  pointwise (12:45)  
```
The above spellings are correct. (The first three are names,
the others I have checked on www.leo.org.)
```

* Found the following (possibly) invalid URLs:  
  URL: http://dx.doi.org/10.1007/978-3-319-50986-0  
    From: man/kader.Rd  
    Status: 404  
    Message: Not Found
```
The error message is presently correct, because the
publication to which the DOI-URL refers is going to
appear not until Oct. 2017.
```

* Non-standard file/directory found at top level:  
  'cran-comments'
```
The directory contain this file.
```



## Downstream dependencies
There are currently no downstream dependencies for this
package.
