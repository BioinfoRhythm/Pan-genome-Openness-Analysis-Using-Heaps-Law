# Pan-genome Openness Analysis Using Heaps' Law

## Overview

This Python script estimates pan-genome openness by modelling cumulative gene-family accumulation using Heaps' law. It uses a gene presence–absence matrix to evaluate pan-genome expansion across randomized genome orders.

The script fits a power-law model to each genome-order iteration and calculates parameter estimates and goodness-of-fit statistics. It also generates a plot showing pan-genome accumulation and core-genome trends.

## Methodology

The pan-genome accumulation curve is modelled using:

\[
P(N)=\kappa N^\gamma
\]

where:

- \(P(N)\) is the cumulative number of gene families identified after analysing \(N\) genomes.
- \(N\) is the number of genomes analysed.
- \(\kappa\) is the scaling parameter.
- \(\gamma\) is the power-law exponent used to characterize pan-genome expansion.

The script performs random genome-order iterations to account for variation associated with genome sampling order. For each iteration, it calculates cumulative gene-family counts and fits the power-law model.

The output includes mean parameter estimates, 95% Student's *t*-confidence intervals across genome-order iterations, the coefficient of determination (\(R^2\)), and root mean square error (RMSE).

A positive estimated \(\gamma\) is interpreted as evidence of an open pan-genome under this model.

## Requirements

- Python 3
- NumPy
- Matplotlib
- SciPy

Install the required packages using:

```bash
pip install numpy matplotlib scipy
```

## Input

The script requires a tab-separated gene presence–absence matrix, such as a Roary `gene_presence_absence.Rtab` file.

- Rows represent gene families.
- Columns represent genomes.
- Values indicate gene presence (1) or absence (0).
- The first column contains gene-family identifiers.

## Usage

```bash
python pan_genome_analysis_corrected.py \
    --input_file gene_presence_absence.Rtab \
    --iterations 10 \
    --output_dir results \
    --threads 30
```

Arguments:

| Argument | Description |
|---|---|
| `-f`, `--input_file` | Input gene presence–absence matrix |
| `-i`, `--iterations` | Number of random genome-order iterations |
| `-o`, `--output_dir` | Directory for output files |
| `--threads` | Thread-count argument; currently not used for parallel computation |

## Output

The output directory contains:

- `pan_core_genome_plot.png` — pan-genome and core-genome accumulation plot.
- `pan_core_genome_plot.pdf` — publication-quality PDF plot.
- `pan_genome_statistics.txt` — parameter estimates, goodness-of-fit statistics, and 95% confidence intervals.

## Citation

If you use this script, please cite the relevant methodological references for Heaps' law and the associated pan-genome analysis approach:

Please add the complete bibliographic details and the original software or repository attribution where applicable.
