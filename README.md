# Code for the UC model with SCC and DR constraints

**PLEASE NOTE, the main models and methodologies are in the listed papers here. Fully understanding these works is the foundation of our work.**
- Short-Circuit Current (SCC) models refer to:
  1. Chu, Zhongda, and Fei Teng. ["Short circuit current constrained UC in high IBG-penetrated power systems." IEEE Transactions on Power Systems 36.4 (2021): 3776-3785.](https://ieeexplore.ieee.org/abstract/document/9329077)
  2. Chu, Zhongda, Jingyi Wu, and Fei Teng. ["Pricing of short circuit current in high IBR-penetrated system." Electric Power Systems Research 235 (2024): 110690.](https://www.sciencedirect.com/science/article/pii/S0378779624005765)
- Data used in this work and relevant work/data please refer to our previous work:
  1. Wang, Peng, and Luis Badesa. ["Imperfect Competition in Markets for Short-Circuit Current Services." arXiv preprint arXiv:2508.09425 (2025)](https://arxiv.org/pdf/2508.09425).
  2. Wang, Peng, and Luis Badesa. ["Pricing Short-Circuit Current via a Primal-Dual Formulation for Preserving Integrality Constraints." arXiv preprint arXiv:2510.05293 (2025)](https://arxiv.org/pdf/2510.05293).
     


**GUIDANCE abot how to use the code of our work**

The work is mainly made of two parts:
1. Modelling SCC constraints.
2. Modelling of DR constraints.

We try to guide you to understand our logistics of coding, once you fully understand, then analyze any power systems you want.
- For the code of SCC modelling, please refer to the files named "_admittance_matrix_calculation.jl_", "_dataset_gene.jl_" and "_offline_trainning.jl_".

  1. "_admittance_matrix_calculation.jl_" calculates the impedance of transmission lines of the system, easy to follow.

  2. "_dataset_gene.jl_" generates the data for classification, i.e., the offline trainning process. The subfunction "_admittance_matrix_calculation.jl_" is called here to obtain the transmission line admittance matrix which is combined with the generators' admittance matrix. In "_dataset_gene.jl_", code from line 79-86 is the equation of actual, exact SCC representation. The remainder of code is generating all possible UC status pairs of generators. The matrix "I_SCC_all_buses_scenarios" are the SCC corresponding to "matrix_ω" (storing UC status and capacity factor of IBR, comprising all possible scenarios).

  3. "_offline_trainning.jl_" is the trainning process, with inputting parameters from above subfunctions.

- For the code of DR constraints, please refer to the file named "_SCC_DR.jl_". It is in the line 252-299, which is easy to follow.
----

If you find something helpful or use this code for your own work, please cite this paper:
<ol>
    Wang, Peng, Zhengmao Li, and Luis Badesa. "Analyzing the Impact of Demand Response on Short-Circuit Current via a Unit Commitment Model." arXiv preprint arXiv:2511.00296 (2025).
</ol>
      <br>
      
<ol> 

@article{wang2025analyzing,<br>
  title={Analyzing the Impact of Demand Response on Short-Circuit Current via a Unit Commitment Model},<br>
  author={Wang, Peng and Li, Zhengmao and Badesa, Luis},<br>
  journal={arXiv preprint arXiv:2511.00296},<br>
  year={2025} <br>
}
