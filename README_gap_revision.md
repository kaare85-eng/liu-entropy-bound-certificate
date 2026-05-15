# Liu Entropy Bound Certificate

This repository contains a draft manuscript and computational certificates for the paper:

**“A Certified Bivariate Entropy Inequality Related to Liu's Entropy Bound”**
by Kåre.

## Current status

This repository is under revision after feedback from Jingbo Liu.

The computational certificates verify a bivariate entropy inequality related to Liu's conditionally i.i.d. coupling method for the union-closed sets conjecture. An earlier draft claimed that this inequality implied the constant

```text
c' = 0.382709087918735...
```

for Frankl's union-closed sets conjecture. That implication is not currently justified.

The issue is distributional: in Liu's conditionally i.i.d. framework, the independent entropy term and the conditionally independent entropy term are averaged with respect to different measures. The pointwise bivariate inequality certified here controls both terms on the same pair `(x, y)`. Therefore, an additional distributional argument or a different inequality is needed to obtain the full union-closed consequence.

## What is certified

The repository currently certifies:

- a bivariate entropy inequality suggested by Liu's optimizer structure;
- a diagonal Bernstein/Rolle certificate;
- off-diagonal Taylor-model determinant certificates;
- branch isolation for the critical-point analysis of the bivariate inequality.

It does **not** currently certify the union-closed lower bound.

## Manuscript

The compiled manuscript is available at:

```text
paper/manuscript.pdf
```

The LaTeX source is available at:

```text
paper/manuscript.tex
```

## Repository structure

```text
liu-entropy-bound-certificate/
├── README.md
├── LICENSE
├── requirements.txt
├── paper/
│   ├── manuscript.pdf
│   └── manuscript.tex
├── code/
│   └── run_all_certificates.py
└── logs/
    └── certificate_output.txt
```

## Reproducing the certificates

Install the only dependency:

```bash
pip install mpmath
```

Run the certificate suite:

```bash
python code/run_all_certificates.py
```

Expected output, last lines:

```text
A.6 Branch isolation: CERTIFIED
---
ALL CERTIFICATES PASSED
Total runtime: ~25 seconds
```

The certificates correspond to Propositions A.1-A.9 in the manuscript. They certify the bivariate inequality and its critical-point analysis, not the union-closed corollary.

## Certificate summary

| Proposition | Content | Method |
|------------|---------|--------|
| A.1 | Constant enclosures | Interval Newton |
| A.2 | Bernstein positivity, min coefficient > 4.55 | Bernstein basis |
| A.3 | Diagonal interval signs | Interval evaluation |
| A.4 | Unique bifurcation point x0 | Interval Newton 1D |
| A.5 | Fold uniqueness, Gamma_u > 0.41 | Interval Newton 2D |
| A.6 | Branch isolation for the bivariate critical-point problem | Topological + band certificate |
| A.7 | Lower arm det D^2 Gamma < 0 | Taylor model |
| A.8 | Upper arm det D^2 Gamma < 0 | Taylor model in (epsilon, s) |
| A.9 | Endpoint Gamma_u > 0 | Analytic bound |

## References

- Frankl (1979) — original conjecture
- Gilmer (2022), arXiv:2211.09055 — entropy method breakthrough
- Liu (2023), arXiv:2306.08824 — conditionally i.i.d. coupling

## Repository

https://github.com/kaare85-eng/liu-entropy-bound-certificate

## Note on development

This paper was developed with substantial AI assistance through extended interaction with large language models. The mathematical arguments and certificates can be verified independently of how they were found.

## License

MIT License. See `LICENSE`.
