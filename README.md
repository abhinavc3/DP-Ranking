# DP-Ranking

Code to reproduce simulation results, figures, and real data analysis results from the paper "Optimal Differentially Private Ranking from Pairwise Comparisons" 

Description of File in the Edge-DP folder
---

## 📂 Main Scripts

- **Simulation-Study-1.R**  
  First baseline simulation study to test Edge-DP under simple settings.

- **Simulation-Study-2.R**  
  Extended study with more parameters varied (dimension, noise level).  
  **Simulation-Study-2 - parallelized.R** provides a parallelized implementation for efficiency on multi-core machines.

- **Simulation-Study-3.R**, **Simulation-Study-4.R**, **Simulation-Study-5.R**  
  Further simulation studies, each exploring a new dimension of the problem (e.g., different graph structures, sampling regimes, or privacy parameters).

- **Experiment-1,2,3.R**  
  Unified script to reproduce the three main experimental setups described in the accompanying paper/manuscript.

- **Helper-Script.R**  
  Provides utility functions (data generation, Edge-DP mechanism implementation, error metrics) used by the other simulation scripts.

- **Plotter.R**  
  Collects raw output from simulations and generates figures for error comparisons (Hamming, L₂, L∞) across settings.

---

## 📂 HPC (Cluster) Code

Located in `hpcc code/`:
- **Experiment-1-2-3-with-L-2(hpcc version).R**  
  HPC-adapted version of the main experiments, optimized for large-scale runs with L₂ error included.
- **log.txt** → Cluster job logs.
- **hpcc-simulation-results/** → Contains large-scale results, organized by study:
  - `accuracy-vs-epsilon/`
  - `accuracy-vs-n/`
  - `accuracy-vs-p/`

Each subfolder includes `.csv` files with simulation outputs for **Hamming error**, **L₂ error**, and **L∞ error** (both estimated and counted).

---

## 📂 Plots

Located in `Plots/`:
- High-level plots of loss functions:
  - **l infinity loss.pdf**
  - **l 2 loss.pdf**
  - **hamming loss via estimation.pdf**
  - **hamming loss via counting.pdf**

- Organized subfolders:
  - **Plots that look good/** → Curated versions of plots, plus supporting parameter file (`theta.txt`).
  - **accuracy-vs-epsilon/**, **accuracy-vs-n/** → Breakdown of how performance scales with privacy budget (ε) and sample size (n).  
    Includes parameter-specific subfolders (e.g., `p=0.5 n=300`, `p=1 epsilon=0.5`).

---


## 🔑 Usage Notes

1. Run `Simulation-Study-*.R` scripts locally for small-scale tests.  
   Use `Simulation-Study-2 - parallelized.R` if working on a machine with multiple cores.  
   For HPC usage, switch to `hpcc code/`.

2. After simulation, run `Plotter.R` to generate figures.  

3. Browse results in `Plots/` (ready-made PDFs) and `hpcc-simulation-results/` (raw `.csv` data).

4. Ignore `.Rhistory` and `.DS_Store` — they are not relevant for analysis.



# Individual-DP Folder

This folder contains code and notebooks for experiments and simulations related to **Individual Differential Privacy (Individual-DP)**.

---

## 📂 Main Files

- **pairwise_comparison.py**  
  Python script implementing pairwise comparison mechanisms under Individual-DP.  
  Contains core functions to run the DP algorithm and evaluate performance.

- **individual_simulate.ipynb**  
  Jupyter notebook for running **synthetic simulations** under Individual-DP.  
  Demonstrates how different parameter settings (e.g., ε, n, noise level) affect outcomes.

- **individual_real_data.ipynb**  
  Jupyter notebook for applying Individual-DP methods to **real datasets**.  
  Useful for validating the methodology beyond synthetic simulations.

---


## 🔑 Usage Notes

1. **Start with the Notebooks**  
   - Use `individual_simulate.ipynb` to explore simulation results and sanity check the methodology.  
   - Then move to `individual_real_data.ipynb` to replicate experiments on real data.

2. **Core Functions**  
   - The `pairwise_comparison.py` script holds reusable functions that can be imported into the notebooks.

