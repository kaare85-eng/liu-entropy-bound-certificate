# Liu Entropy Bound Certificate

This repository contains the manuscript and computational certificates
for the paper:

**“A Direct Proof of Liu’s Entropy Bound for the Union-Closed Sets Conjecture”**
by Kåre

## What the paper proves

We prove a bivariate entropy inequality which rigorously establishes
the constant

```
c' = 0.382709087918735...
```

for Frankl’s union-closed sets conjecture. This constant was obtained
conditionally by Liu (2023, arXiv:2306.08824). Our proof replaces
Liu’s numerical optimizer-structure hypothesis with a direct certified
proof of the underlying entropy inequality.

## Repository structure

```
liu-entropy-bound-certificate/
├── README.md
├── LICENSE
├── requirements.txt
├── paper/
│   ├── manuscript.pdf        (compiled paper)
│   └── manuscript.tex        (LaTeX source)
├── code/
│   └── run_all_certificates.py
└── logs/
    └── certificate_output.txt
```

## Reproducing the certificates

Install the only dependency:

```
pip install mpmath
```

Run the certificate suite:

```
python code/run_all_certificates.py
```

Expected output (last lines):

```
A.6 Branch isolation: CERTIFIED
---
ALL CERTIFICATES PASSED
Total runtime: ~25 seconds
```

The certificates correspond to Propositions A.1–A.9 in the paper.

## Certificate summary

|Proposition|Content                          |Method                 |
|-----------|---------------------------------|-----------------------|
|A.1        |Constant enclosures              |Interval Newton        |
|A.2        |Bernstein positivity (min > 4.55)|Bernstein basis        |
|A.3        |Diagonal interval signs          |Interval evaluation    |
|A.4        |Unique bifurcation point x₀      |Interval Newton 1D     |
|A.5        |Fold uniqueness, Γ_u > 0.41      |Interval Newton 2D     |
|A.6        |Branch isolation                 |Topological + band cert|
|A.7        |Lower arm det D²Γ < 0            |Taylor-model           |
|A.8        |Upper arm det D²Γ < 0            |Taylor-model in (ε,s)  |
|A.9        |Endpoint Γ_u > 0                 |Analytic bound         |

## References

- Frankl (1979) — original conjecture
- Gilmer (2022) arXiv:2211.09055 — entropy method breakthrough
- Liu (2023) arXiv:2306.08824 — conditionally i.i.d. coupling

## Repository
https://github.com/kaare85-eng/liu-entropy-bound-certificate

## Note on development

This paper was developed with substantial AI assistance through
extended interaction with large language models. The mathematical
arguments and certificates can be verified independently of how
they were found.

## License

MIT License. See LICENSE file.
