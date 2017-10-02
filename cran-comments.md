---
title: "cran-comments.md"
author: "Gerrit Eichner"
date: "October 2, 2017"
output:
  pdf_document: default
  html_document: default
---

## Resubmission
This is a resubmission. In this version I have (as requested):  
* added at least one executable example in each .Rd-file.

## Test environments
* local: Windows 10 Pro on x86_64-w64-mingw32/x64 (64-bit), R 3.4.2
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.4.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE (from win-builder):

* Possibly mis-spelled words in DESCRIPTION:  
  Eichner (12:48, 14:32)
  Srihera (11:73)
  Stute (11:63, 12:60, 14:44)
  nonparametric (10:67)
  pointwise (13:42)
```
The above spellings are correct. (The first three are
names, the others I have checked on www.leo.org.)
```

## Downstream dependencies
There are currently no downstream dependencies for this
package.
