# FMS

## Preface

This repository implements command line interfaces to compute:
- $p_{n, k}$ for k in [0, w - 1], given parameters n and w
- Pr[Z = k] for k in [0, (w - 1)*m], given parameters n, w and m
- Find the smallest $N_0$ satisfying the differential privacy, given parameters m, w, $\epsilon$ and $\delta$.

## Requirements
- OS: Ubuntu 16.04 LTS
- GCC, version 7.4.0
- [Arb](http://arblib.org/)

## To compile this repository
Download this repository and run the following commands in shell.
```bash
cd fms/
make all
```


## Command line arguments

**Warning: This program is not optimized and may take a long time to compute the result.**

| Parameters                    | Description
| ---------                     | -------- 
| -i                            | Works to compute the PGF of z. -n and -w should be set.
| -a                            | Works to compute the PMF of Z. -n, -w and -m should be set.
| -s                            | Works to find the smallest $N_0$ satisfying the differential privacy. -w, -m, -e and -d should be set. -c flag is optional to speed up the computation.
| -n [n]                        | Value of n.
| -w [w]                        | Value of w.
| -m [m]                        | Value of m.
| -e [epsilon]                  | Value of $\epsilon$ (in float type).
| -d [delta]                    | Value of $\delta$ (in int type). E.g., -d 40 represents $\delta = 2^{-40}$.
| -o [file]                | Specify the output file.
| -c [ceiling]                  | Ceiling of the interval in which our binary search is performed to find the smallest $N_0$. This ceiling will be automatically doubled if there is no such $N_0$ in current interval. Default 8192.

