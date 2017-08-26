---
title: "cran-comments.md"
author: "Gerrit Eichner"
date: "26 August 2017"
output: html_document
---

## Resubmission
This is a resubmission. In this version I have (as requested):  
* added DOIs for the references in the DESCRIPTION.  
* shortened or "switched off" (\donttest{}) the examples which
  were running > 5 sec. (Now all run in < 5 sec.) 

## Test environments
* local: Windows 10 Pro on x86_64-w64-mingw32/x64 (64-bit), R 3.4.1
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.4.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE (from win-builder):

* Possibly mis-spelled words in DESCRIPTION:  
  Eichner (12:16, 13:32)  
  Srihera (11:73)  
  Stute (11:63, 12:28, 13:44)  
  nonparametric (10:67)  
  pointwise (12:45)  
```
The above spellings are correct. (The first three are
names, the others I have checked on www.leo.org.)
```

## Downstream dependencies
There are currently no downstream dependencies for this
package.
