<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/Quantile_Coherence_Clustering">
    <img src="UoY.png" alt="York Logo" height="150">
    <img src="Kaustlogo.png" alt="KAUST Logo" height="150">
    <img src="IBM.png" alt="IBM Logo" height="150">
  </a>

<h3 align="center">A semi-parametric estimation method for quantile coherence with an application to bivariate financial time series clustering</h3>
</div>

## Abstract
<p align="justify">
In multivariate time series analysis, spectral coherence traditionally measures linear dependencies between time series across different frequencies but often fails to capture nonlinear dependencies. In contrast, quantile coherence detects these nonlinear relationships across various quantile levels using trigonometric quantile regression. A new semi-parametric technique for estimating quantile coherence is introduced, combining the parametric spectrum of a vector autoregressive (VAR) model with nonparametric smoothing across quantiles. For each quantile level, the quantile autocovariance function (QACF) is obtained by applying the Fourier inverse transform to quantile periodograms. The multivariate Durbin-Levinson algorithm is then used to estimate the VAR parameters, which are subsequently applied to derive the quantile coherence estimate. A nonparametric smoother is applied across quantiles to enhance the initial estimate. Numerical results demonstrate that this method outperforms traditional nonparametric approaches. Moreover, clustering bivariate time series based on quantile coherence shows advantages over using ordinary VAR coherence. For instance, when applied to financial stocks, quantile coherence identifies clusters with a more informative structure, providing insights into diversified investment portfolios, which could help investors make better decisions.
</p>

### Main Results
The R script files in the `R Code` folder should be used in the following order:

#### 1. Setup, Package Initialization, and Auxiliary Utilities
* **aux_semi_parametric_estimation_QC.R**: Houses foundational backend auxiliary functions required for computing the semiparametric estimator of quantile coherence.
* **fviz_nbclust.R**: A customized version of the `fviz_nbclust` utility from the `NbClust` library, optimized for standardized axis scales and clean label formatting.
* **Financial_data.csv**: The main empirical dataset housing the raw bivariate financial time series used throughout the paper's stock market case study.
* **beta.rds** & **stocks_names.rds**: Serialized R data storage files holding pre-calculated model parameters ($\beta$ CAPM) and structural financial stock metadata.

#### 2. Quantile Coherence Semi-Parametric Estimation Frameworks
* **Semi_parametric_estimator_QC_SD.R**: Implements the primary semiparametric estimator of quantile coherence tailored specifically for simulated data testing.
* **Semi_parametric_estimator_QC_FTS.R**: Deploys the main semiparametric estimator framework optimized for real-world empirical Financial Time Series.
* **Joint_smoothing.R**: Implements the nonparametric joint smoothing splines simultaneously across multiple quantile levels to stabilize the initial matrix profiles.

#### 3. Baseline Validation, Benchmarks, and Competing Models
* **VAR_Coh.R**: Computes the theoretical/parametric ordinary Vector Autoregressive coherence metrics used as linear model benchmarks.
* **Alternative_methods.R**: Executes alternative or competing legacy estimation and clustering baseline regimes to test against the proposed semiparametric strategy.

#### 4. Financial Clustering Framework and Applications
* **Clustering_simulation.R**: Runs controlled simulation pipelines evaluating the accuracy, stability, and group assignments of the proposed clustering technique.
* **Computations_OC.R**: Generates benchmark financial stock partitions based entirely on ordinary linear VAR coherence profiles.
* **Clustering_FTS_QC.R**: Runs the master functional time series clustering algorithm on the empirical stock dataset using nonlinear quantile coherence matrices as the core distance/dissimilarity metric.

## How to Cite
If you use this code, data, or methodology in your research, please cite our paper:

**APA Style:**
> Jiménez-Varón, C. F., Sun, Y., & Li, T.-H. (2024). A semi-parametric estimation method for quantile coherence with an application to bivariate financial time series clustering. *Econometrics and Statistics*. https://doi.org/10.1016/j.ecosta.2024.11.002

**BibTeX:**
```bibtex
@article{jimenezvaron2024semiparametric,
  author    = {Jim{\'e}nez-Var{\'o}n, Cristian F. and Sun, Ying and Li, Ta-Hsin},
  title     = {A semi-parametric estimation method for quantile coherence with an application to bivariate financial time series clustering},
  journal   = {Econometrics and Statistics},
  year      = {2024},
  issn      = {2452-3062},
  doi       = {10.1016/j.ecosta.2024.11.002},
  url       = {[https://doi.org/10.1016/j.ecosta.2024.11.002](https://doi.org/10.1016/j.ecosta.2024.11.002)}
}
