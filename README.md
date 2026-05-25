# dwell-time-summing

Supporting code and data for:

> **Processing Ion Currents using Multiple Dwell Times**
>
> Fred Lytle<sup>1</sup>, Sangeeta Pandey<sup>2</sup>, Randall K. Julian, Jr.<sup>1</sup>\*
>
> 1. Indigo BioAutomation, Carmel, Indiana, United States
> 2. Department of Chemistry, Purdue University, West Lafayette, Indiana, United States
>
> \*Email: rkjulian@indigobio.com

## Overview

This repository contains the raw data and R analysis scripts for a study investigating whether averaging rates over individual dwell times in LC-MS/MS improves concentration precision relative to a single measurement over the largest possible dwell time. The chemical system is naloxone/naloxone-d5 measured on a Waters Xevo TQD at five concentrations in both flow injection and chromatographic modes, using three timing configurations: a single 96 ms dwell, four 21 ms dwells, and ten 6 ms dwells (each summing to a 100 ms cycle time with 4 ms pauses).

Rate averaging and count summation give statistically indistinguishable coefficients of variation, and neither improves concentration precision relative to the single long dwell. The chromatographic signal-to-noise ratio (peak height / baseline noise) increases under summation, but this statistic measures peak detectability, not concentration precision. True signal-to-noise (1/CV) is unchanged. The reported LLOQ gains in the literature arise from a detectability criterion embedded in dual-criterion LLOQ definitions (S/N > 10 and CV < 20%), not from improved quantitative precision.

## Repository Structure

```
code/                        # R analysis scripts (snake_case)
  functions/                 # Shared R functions (EMG fitting, optimization)
data/
  flow_injection_data/       # Flow injection raw data (2021-08-23)
  chromatogram_data/         # Chromatography raw data (2021-08-25)
```

## Running the Analysis

All R scripts assume the working directory is the repository root. Scripts read data from `data/` using structured filenames to select by concentration, dwell time, replicate, and transition.

### Requirements

- R (>= 4.0)
- CRAN packages: `MASS`, `NADA`
- No Bioconductor or tidyverse dependencies

### Script-to-Manuscript Map

| Manuscript Element | Script(s) |
|---|---|
| Figure 1 (Timing diagram) | `timing_figure.R` |
| Figure 2 (Flow injection traces) | `flow_injection_full_trace_rate.R`, `flow_injection_full_trace_count.R` |
| Figure 3 (Excess noise CV vs concentration) | `dwell_time_cv_vs_concentration.R` |
| Figure 4 (Combined calibration, panels A+B) | `flow_injection_calibration_combined.R` |
| Figure 5 (Chromatogram peak) | `chromatogram_peak_counts.R`, `chromatogram_peak_rates.R` |
| Table 1 (Plateau statistics) | `flow_injection_full_trace_rate.R`, `flow_injection_full_trace_count.R`, `table_1_values.R` |
| Table 2 (Flow injection calibration) | `flow_injection_calibration_rate.R`, `flow_injection_calibration_count.R`, `table_2_slope_corrections.R` |
| Table 3 (EMG peak areas) | `chromatogram_peak_rates.R`, `chromatogram_peak_counts.R` |
| Table 4 (Chromatography calibration) | `chromatogram_calibration_peak_rates.R`, `chromatogram_calibration_peak_counts.R` |
| Table 5 (Correlation and covariance) | `flow_injection_dt_rate_average_correlations.R`, `flow_injection_dt_count_average_correlations.R` |
| Table S1 (Reverse calibration, flow injection) | `flow_injection_reverse_calibration_rate_simulation.R`, `flow_injection_reverse_calibration_count_simulation.R` |
| Table S2 (Reverse calibration, chromatography) | `chromatogram_reverse_calibration_rate_simulation.R`, `chromatogram_reverse_calibration_count_simulation.R` |
| Equation 4 (Covariance demonstration) | `data_with_covariance_computations.R`, `variance_of_an_average_of_correlated_variables.R` |
| CV confidence intervals | `cv_confidence_intervals.R` |
| Dwell time corrections | `dwell_time_corrections.R`, `table_2_slope_corrections.R` |
| Ion-arrival simulation | `dwell_time_simulation_counts.R`, `dwell_time_simulation_rates.R` |
| Baseline noise characterization | `baseline_detector_noise.R` |
| Excess noise vs bandwidth | `naloxone_excess_noise_versus_bandwidth.R` |
| IS proportional noise | `naloxone_d5_proportional_noise_cv_vs_dwell_time.R` |
| 96 ms autocorrelation | `flow_injection_dt96_correlation_investigation.R` |
| Raw vs detrended correlations | `table7_raw_vs_detrended_comparison.R` |
| Plateau drift visualization | `flow_injection_plateau_slope_histograms.R` |

### Shared Functions (`code/functions/`)

| File | Purpose |
|------|---------|
| `EMG.R` | Exponentially modified Gaussian for chromatographic peak fitting |
| `marquardtSearch.R` | Levenberg-Marquardt nonlinear least-squares optimizer |
| `getChisqr.R` | Chi-square objective function |
| `getDeriv.R` | Numerical derivatives for the optimizer |
| `getIndexSelectedFilenameFlow.R` | Parse structured filenames (flow injection data) |
| `getIndexSelectedFilenameChrom.R` | Parse structured filenames (chromatogram data) |

## Raw Data

### Flow Injection (2021-08-23)

Two-column text files (time in seconds, rate in counts/sec). Three replicates at each of five concentrations (50, 250, 500, 1000, 1500 ng/mL) for naloxone and naloxone-d5, acquired with three dwell time configurations (6, 21, 96 ms) and 4 ms pause time. Filename encodes all acquisition parameters:

```
{date}_{conc}_Naloxone+{IS}_{channels}_DT_{dwell}ms_PT_{pause}ms_{rep}_{transition}_{energy}_Pos.txt
```

### Chromatography (2021-08-25)

Same format and naming convention. Four replicates at five concentrations (50, 100, 200, 250, 500 ng/mL). EMG peak fitting applied by the chromatogram analysis scripts.

## Experimental Parameters

- **Analyte:** naloxone, m/z 328.1 -> 212.1, collision energy 40 eV
- **Internal standard:** naloxone-d5, m/z 333.1 -> 212.1, collision energy 36 eV
- **Instrument:** Waters Xevo TQD
- **Timing configurations:** 1x96 ms, 4x21 ms, 10x6 ms (cycle time ~100 ms, pause time 4 ms)
- **Flow injection concentrations:** 50, 250, 500, 1000, 1500 ng/mL (3 replicates)
- **Chromatography concentrations:** 50, 100, 200, 250, 500 ng/mL (4 replicates)

## License

Apache License 2.0. See [LICENSE](LICENSE) for details.
