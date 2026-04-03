# Zeus_Notebook

Public notebook: **[`notebooks/gcft_hera_f2_xion_prior.ipynb`](notebooks/gcft_hera_f2_xion_prior.ipynb)**

**Data**

- [HEPData ins836107](https://www.hepdata.net/record/ins836107) — H1+ZEUS **2010** combined inclusive NC $e^+p$ ($F_2$, $\sigma_r$), Aaron et al., JHEP 01 (2010) 109.
- [HEPData ins1377206](https://www.hepdata.net/record/ins1377206) — **2015** HERA combined NC $\sigma_r$ at $\sqrt s=318$ GeV (secondary coverage plot in the notebook; not merged with 2010 $F_2$).

**Notebook features**

- Merged $F_2$ + $\sigma_r$; **stat** vs **quadrature-combined** relative uncertainties from all `%` columns on each $\sigma_r$ row (diagonal / independent-component approximation — **not** the full YAML covariance).
- **Global** log-plane and **global** quadratic surface in $(\log x,\log Q^2)$.
- **Per-$Q^2$** fit: **quadratic** in $\log x$ if $\geq 4$ points, else **line**, else global fallback.
- Illustrative **xion** prior samples; **MAP** on amplitude $\eta$ (scipy) at fixed $x_0$ with diagonal $\chi^2$ + Gaussian prior.
- Residuals vs $Q^2$, **dual** pull histogram (comb vs stat $\sigma$), $\sigma_r(x)$ demo, 2010 and 2015 **$(Q^2,x)$** coverage plots, four **$F_2(x)$** panels (`pick_q2_slices` by point count).
- Links to **YAML** bundles for ins836107 / ins1377206; `pyyaml` is installed/imported in Colab for follow-on parsing.

**Colab:** [open in Colab](https://colab.research.google.com/github/ChaonickGCFT/Zeus_Notebook/blob/main/notebooks/gcft_hera_f2_xion_prior.ipynb).

**Not included:** dedicated $F_L$ tables, diffractive channels, CC — add other HepData records and a YAML→covariance step for a publication-grade fit.
