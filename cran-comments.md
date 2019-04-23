#  Version 0.1.0

## Test environments
* local OS X install, R 3.5.0
* ubuntu 14.04 (on travis-ci), R 3.5.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 2 notes

* One note refers to mis-spelled words. 

```
checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Maikol Solís <maikol.solis@ucr.ac.cr>’
```

```
Possibly mis-spelled words in DESCRIPTION:
Sobol (3:22, 7:40)
Solís (7:217)
```

* Second note refers to time to execution in the examples (~ 2.5 minutes)
```
Examples with CPU or elapsed time > 5s
        
         user system elapsed
plot    46.66   0.05   46.77
print   46.55   0.05   46.87
sobolnp 45.72   0.08   45.87
```
