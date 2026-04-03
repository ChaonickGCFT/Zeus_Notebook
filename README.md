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
- **v0.3:** fixed-$Q^2$ $F_2$ **residual slice** vs $\log_{10} x$ (table $Q^2$ nearest **4** GeV$^2$, tie → **4.5** GeV$^2$) and **`hera_systematics_quadrature.csv`** — stat, uncorrelated experimental, correlated experimental, three combination-correlated sources, and overall normalisation as **separate relative columns** (from `sys_1`–`sys_5` + norm on each $\sigma_r$ row per the [ins836107](https://www.hepdata.net/record/ins836107) record text).

**Colab:** [open in Colab](https://colab.research.google.com/github/ChaonickGCFT/Zeus_Notebook/blob/main/notebooks/gcft_hera_f2_xion_prior.ipynb).

**Changelog**

- **v0.3** — `f2_residual_slice_z.png` + `hera_systematics_quadrature.csv` under [`notebooks/v0.3_assets/`](notebooks/v0.3_assets/); notebook cells with raw CSV link; loader in `notebooks/_v02_generate_outputs.py` extended to keep $\sigma_r$ systematics split for merged $F_2$ rows. Raw systematics: `https://raw.githubusercontent.com/ChaonickGCFT/Zeus_Notebook/main/notebooks/v0.3_assets/hera_systematics_quadrature.csv`
- **v0.2** — Freeze table $Q^2 > 12k$ → $12k$ template; residual heat-map; `chi2_per_q2.csv`; diffractive template CSV; short MCMC trace.

**v0.2 (freeze-$Q^2$ + release assets)**  
The top of the $Q^2$ ladder is **stabilised** by copying the per-$Q^2$ $(\log x)$ template from the largest table $Q^2 \le 12\,000\,\mathrm{GeV}^2$ onto any table with **higher** $Q^2$ (for ins836107 that is only the $Q^2 = 20\,000\,\mathrm{GeV}^2$ slice — lower table $Q^2$ keep their own fits). The residual heat-map uses a **$\pm 7\%$** colour scale; most of the $(\log Q^2,\log x)$ plane stays inside that band, with the **highest-$Q^2$ column** showing the largest mean $|\Delta F_2|/F_2$. **$\Sigma\chi^2$** (diagonal combined $\sigma$) is **largest** for the medium-$Q^2$ tables ($\sim 500$–$1500\,\mathrm{GeV}^2$); see `chi2_per_q2.csv`. Global **median** $|\Delta F_2|/F_2$ remains $\mathcal{O}(2\%)$ with $\sim 77\%$ of points below $7\%$ (regenerate via `python notebooks/_v02_generate_outputs.py` from repo root).

Committed outputs under [`notebooks/v0.2_assets/`](notebooks/v0.2_assets/):

- `residual_heatmap_v02.png` — mean relative residual in $(\log_{10} Q^2,\log_{10} x)$ (baseline, no xion)  
- `chi2_per_q2.csv` — $\sum\chi^2$, $n$, mean $|\Delta F_2|/F_2$, and fraction below $7\%$ **per table $Q^2$**  
- `chi2_top_bins.csv` — $\sum\chi^2$ by wide $Q^2$ bin (legacy summary)  
- `sigma_diff_gcft_vs_data.csv` — H1 diffractive-style $x F_2^{D(3)}$ from [HEPData ins447269](https://www.hepdata.net/record/ins447269) plus a **log–log template** column `model_gcft_template` (phenomenological cross-check, not a full GCFT fit). Raw URL for tools: `https://raw.githubusercontent.com/ChaonickGCFT/Zeus_Notebook/main/notebooks/v0.2_assets/sigma_diff_gcft_vs_data.csv`  
- `mcmc_eta_trace.npz`, `mcmc_eta_trace.png` — short Metropolis trace for the illustrative xion amplitude $\eta$

**v0.3 assets** — [`notebooks/v0.3_assets/`](notebooks/v0.3_assets/)

- `f2_residual_slice_z.png` — baseline $F_2$ residuals vs $\log_{10} x$ at $Q^2 \approx 4.5\,\mathrm{GeV}^2$ (nearest table bin to $4\,\mathrm{GeV}^2$; tie-break upward)  
- `hera_systematics_quadrature.csv` — per-point $\sigma$ breakdown + $z=Q^2/(E_{\rm cm}^2 x)$ for cross-checks with diffractive / DIS kinematics  

**Not included:** dedicated $F_L$ tables, full diffractive theory fit, CC — add other HepData records and a YAML→covariance step for a publication-grade fit.
