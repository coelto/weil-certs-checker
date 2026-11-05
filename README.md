# Independent verification of certificates of lower bounds for the truncated Weil form

Vadim Chekmachev · chekmachev@phystech.edu · ORCID [0009-0009-2982-0532](https://orcid.org/0009-0009-2982-0532)

Paper: [![10.5281/zenodo.23848596](https://zenodo.org/badge/DOI/10.5281/zenodo.23848596.svg)](https://doi.org/10.5281/zenodo.23848596) · Certificate data: [![10.5281/zenodo.23091912](https://zenodo.org/badge/DOI/10.5281/zenodo.23091912.svg)](https://doi.org/10.5281/zenodo.23091912)

A small Python program that re-verifies the certificates from the paper on rigorous lower bounds for the truncated Weil form using the published certificate data. It takes nothing from the program that constructed the certificates except the published symbol balls and the certificate parameters. All computations are carried out in the ball arithmetic of Arb from the FLINT library.

## What is verified

The certificate yields the inequality `Q_c(f) >= eps ||f||^2` in its sector provided that the three conditions of Theorem 5.10 of the paper hold. The program verifies each of them.

Symbol balls. At each retained node the symbol is recomputed from the definition, and its ball must lie inside the published one (`check_symbol.py`).

Threshold outside the retained nodes. The retained nodes must cover the initial segments, the symbol at every other node up to the end of the scan must be no less than the threshold, and the monotone lower bound after it must be no less than the threshold (`check_scan.py`, Lemmas 4.2 and 5.8).

Matrix inequality. The matrix `Id - W_up` is reassembled from the published balls via Proposition 5.2, Lemma 5.3 and Proposition 5.9, and its positive definiteness is proved by a ball Cholesky decomposition (`check_matrix.py`, Lemma 5.4).

## Installation

Requires Python 3.10 or newer and the python-flint package, which is installed together with the FLINT library.

```
pip install -r requirements.txt
```

## Data

The certificate data reside in the Zenodo record indicated below. Both archives should be unpacked into a single folder. This yields a folder `truncated-weil-form-certificates-en` containing `certificates/` and `balls/`. This folder is passed to the program as the first argument.

Certificate data: [![10.5281/zenodo.23091912](https://zenodo.org/badge/DOI/10.5281/zenodo.23091912.svg)](https://doi.org/10.5281/zenodo.23091912)

```
truncated-weil-form-certificates-en-v1.1.zip   md5 5dba08db69804797817a500f9be729fd
truncated-weil-form-balls-en-v1.1.zip          md5 6a4240106d12c81b17712d85d1aff377
```

## Running

A full check of one certificate looks like this.

```
python verify.py truncated-weil-form-certificates-en c10_even_eps
```

The outcome `VERIFIED` means that all three conditions have been verified completely. The outcome `VERIFIED (partial scan)` means that the symbol and the matrix have been verified completely, while the node scan was performed on a sample or on an interval. The outcome `NOT VERIFIED` means that at least one condition fails. The flags are described in the output of `python verify.py --help`.

## Running time and completeness of the verification

The running times below were measured on a single computer with 18 worker processes. A full node scan grows very fast together with `c`. For `c <= 31` the verification is complete and takes from seconds to tens of minutes. For `c = 45`, `47` and `71` a full scan would take days, so the symbol and the matrix are verified completely, while the nodes are checked on a sample or interval by interval, and the intervals can be distributed among computers.

| Certificate | Order of the matrix | Scanned nodes | Time | Outcome |
|---|---:|---:|---:|---|
| `c2_even_eps` | 3 | 0 | 3 s | `VERIFIED` |
| `c5_odd_eps` | 9 | 0 | 3 s | `VERIFIED` |
| `c10_even_eps` | 33 | 9 572 | 4 s | `VERIFIED` |
| `c17_odd_eps` | 187 | 303 459 | 9 s | `VERIFIED` |
| `c20_even_eps` | 311 | 4 996 653 | 40 s | `VERIFIED` |
| `c23_odd_eps` | 457 | 5 219 530 | 45 s | `VERIFIED` |
| `c31_even_eps` | 1663 | 215 197 634 | 33 min | `VERIFIED` |

## Negative controls

The program `controls.py` introduces defects into a temporary copy of the data or the code and checks that each check detects them: a shift of the published ball, doubling of `eps`, a shortened scan, an inflated threshold, an extended initial segment, a wrong sign of the pole constant and doubled atomic weights. Before this, the unmodified certificate must pass all three checks.

## License

The code is distributed under the MIT license, see the file [`LICENSE`](LICENSE). The paper and the data have their own terms, which are stated in their Zenodo records.

## How to cite

Cite the paper by its DOI. The certificate data are cited by the DOI of their record.

* Paper: Rigorous lower bounds for the truncated Weil form up to support length ln 71 via a lattice minorant, [![10.5281/zenodo.23848596](https://zenodo.org/badge/DOI/10.5281/zenodo.23848596.svg)](https://doi.org/10.5281/zenodo.23848596)
* Certificate data: [![10.5281/zenodo.23091912](https://zenodo.org/badge/DOI/10.5281/zenodo.23091912.svg)](https://doi.org/10.5281/zenodo.23091912)
