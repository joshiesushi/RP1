Project title:
Multiscale ODE-CPM Modelling of Cisplatin Toxicity in Proximal Tubule Segments: Incorporating Nephron Cell Heterogeneity

---------------------------------------------------

Author: Carl Joshua Sandico Eugenio

---------------------------------------------------

Introduction:
Despite the effectiveness of cisplatin, its usage is limited, due to nephrotoxicity. It is therefore of great importance to unravel the mechanisms behind Drug-induced Kidney Injury (DIKI). This study employs a three-level multiscale ODE-CPM model to simulate the nephrotoxicity in the nephron, while incorporating the heterogenous nature of the kidney.

---------------------------------------------------

Intallation:
Programs needed to run the codes are the following,

1. Rstudio - https://www.rstudio.com/categories/rstudio-ide/
2. R - https://www.r-project.org
3. cmdstan - https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html
4. Morpheus - https://morpheus.gitlab.io
5. ImageJ - https://imagej.nih.gov

Features:
1. Data Analysis/Final Data Analysis.Rmd
This Rmd file is used to analyse RNA-seq data and to calculate a cellular event module score per replicate. Pre-defined gene sets are used to represent the scores of the modules.
  - Loads RNA-seq data 
  - Performs quality check of the data 
      -Removing low %mapped samples
      -Aggregating technical replicates
      -Removing samples with mean lower than first percentile
      -Removing low correlated replicates per sample
  - Normalizes counts using CMP and transofrming to log2
  - Calculating log2 Fold-change of the data
  - Using predefined list of genes to calculate score estimates for each cellular event module
  - Plotting calculated score
  - Analysis of data


2. Rstan/Cisplatin Exposure
This folder includes a number of R scripts which can be run to perform bayesian paramtere inference using cmdstanR. The fits are then saved automatically with the name stating the number of iterations, the name of the model, and the date of inference.
  - Includes stan models for different types of the ODE model (LDLR, LDMR, MDMR, and MDLR). The variable-dependent production rate of the state variables can either be linear (L) or explained with a Michaelis Menten-like term (M).
  - Each stan model block can be editted to your liking to included constraints and priors as you see fit. 
  -To visualize the fit open either one of the 2 available Rmd files. This shows a step by step how to visualize and run the model with the parameter estimates infereed using cmdstanR.
  - Output is saved in output folder.

3. Morpheus/Final_CPM_Model_MDMR.xml
This XML code can be opened editted. For a user friendly experience I recommend openingg this file with Morpheus. This CPM model depicts the nephrotoxicity in either CPT or PPT cell types of the proximal tubules over time. The initial plasma concentration of platinum can be changed to run in-silico experiments. The current settings with 45 initial plasma value is the conditions for 5 mg/kg cisplatin treatment in rats. Adjust this accordingly.

4. Morpheus Histo/Plotting_Histo.Rmd
This file is used to visualized the data acquired using image analysis using ImageJ.

---------------------------------------------------

Contributions and Acknoledgement:
The completion of this research would not have been possible without the invaluable assistance, guidance, and support provided by my daily supervisor, Filippo di Tillio. The formulation used for the score calculation originated from him. His continuous feedback significantly enriched my research endeavours. Gratitude is also extended to Joost Beltman, whose valuable suggestions and feedback proved useful in accomplishing this study. Furthermore, appreciation is owed to everyone in the Image-based Computational Biology group at LACDR, including Elsje Burgers, Raju Sharma, Levi Winkelman, and Justin Chotoe, for their contributions and input. Lastly, the data employed in this research are attributed to Lukas Wijaya, whose valuable contributions have been an integral part of this project.

---------------------------------------------------

Contact info:
For any question feel free to contact me,
email: eugeniocarl97@gmail.com
