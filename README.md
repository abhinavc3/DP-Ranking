# DP-Ranking

Code to reproduce simulation results, real data analysis results, and figures in the paper "Optimal Differentially Private Ranking from Pairwise Comparisons"

For reproducing experiment results in the paper, you may follow the steps below.

### 4.1 Simulated Data Experiments under Edge DP

Code in `/Edge-DP/script/` produces raw numerical results. These numerical results are then manually collected in four subfolders, each corresponding to one experiment. The notebook `/Edge-DP/edge_plot.ipynb` provides functions for plotting numerical results.

Experiment 1: Run experiment 1 section of `Experiment-1-2-3.R` script, and the results get saved in `Simulation-Results/accuracy-vs-n` folder. Manually transfer them to `/Edge-DP/accuracy-vs-n/*.csv`. Run Section 1 of `/Edge-DP/edge_plot.ipynb` to plot them. Each subsection corresponds to one subfigure of Figure 1.

Experiment 2: Run experiment 2 section of `Experiment-1-2-3.R` script, and the results get saved in `Simulation-Results/accuracy-vs-p` folder. Manually transfer them to  `/Edge-DP/accuracy-vs-p/*.csv`. Run Section 2 of `/Edge-DP/edge_plot.ipynb` to plot them. Each subsection corresponds to one subfigure of Figure 2.

Experiment 3: Run experiment 3 section of `Experiment-1-2-3.R` script, and the results get saved in `Simulation-Results/accuracy-vs-epsilon` folder. Manually transfer them to  `/Edge-DP/accuracy-vs-epsilon/*.csv`. Run Section 3 of `/Edge-DP/edge_plot.ipynb` to plot them. Each subsection corresponds to one subfigure of Figure 3.

Experiment 4: Run `Experiment-4.R` script, the results get saved in `Simulation-Results/one-shot-comparison` folder. Manually convert them to csv and transfer them to `/Edge-DP/one-shot-comparison/*.csv`. Sections 1 and 3 of `/Edge-DP/edge_plot.ipynb` can be modified to produce the two subfigures of Figure 4.


### 4.2 Simulated Data Experiments under Individual DP

The main functions are defined in  `/Individual-DP/pairwise_comparison.py`. The notebook `/Individual-DP/individual_simulate.ipynb.py` loads them and defines additional helper functions for producing experiment results.

Experiment 5: run Section 1 of `/Individual-DP/individual_simulate.ipynb.py`. Each subsection corresponds to one subfigure of Figure 5.
Experiment 6: run Section 2 of `/Individual-DP/individual_simulate.ipynb.py`. Each subsection corresponds to one subfigure of Figure 6.
Experiment 7: run Section 3 of `/Individual-DP/individual_simulate.ipynb.py`. Each subsection corresponds to one subfigure of Figure 7.

### 4.3 Real Data Analysis under Individual DP
The main functions are defined in  `/Individual-DP/pairwise_comparison.py`. The notebook `/Individual-DP/individual_real_data.ipynb` loads them and defines additional helper functions for producing experiment results.

The two subsections of `/Individual-DP/individual_real_data.ipynb` each corresponds to one real data set and produces one subfigure of Figure 8.
