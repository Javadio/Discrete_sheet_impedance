# Discrete_sheet_impedance
---------------------------------------------------------------------------
This repository contains MATLAB code to reproduce some of the results from the paper "Advancing RIS Beamforming Efficiency: Moving
Beyond Diagonal Matrix Techniques and "[Physically Consistent RIS: From Reradiation Mode Optimization to Practical Realization](<https://arxiv.org/abs/2409.17738>). 

# Abstract
---------------------------------------------------------------------------
Optimizing wireless propagation channels is crucial for the advancement of future communication technologies. Reconfigurable Intelligent Surfaces (RIS) — flat panels composed of active or passive elements — offer a transformative approach to reshaping the communication environment.
Currently, the diagonal matrix technique is commonly employed for RIS design in RIS-based wireless communication systems. However, this method faces a significant challenge: low scattering efficiency. To address this issue, we propose a novel approach called the "Discrete Sheet Impedance Model".

This repository includes all the necessary files, enabling users to reproduce the results presented in the referenced papers.
# Results
---------------------------------------------------------------------------
In the file optim8GHz.mat, you can define your frequency, number of the unit cells in the supecell, desired deflection frequency, substrate thickness and all the required information. By defining the cost function as $F({Z_{1...K}}) = \left| {{A_{\rm cal}} - {A_{\rm goal}}} \right|$ and employing the MultiStart and Fmincon optimization algorithms, Matlab searches the proper matrix ${Z_{1...K}}$ that minimizes the cost function i.e. maximal reflected power  in the prescribed direction w.r.t the incident power. Since our aim is to increase the reflected power only for desired deflection angle and minimizing it for all the other harmonics, for an ideal MS, abs(R(N+2)-1/sqrt(cosd(th_d))). 
This Matlab file generates $K$ numbers of ${Z_{1...K}}$. In our specific example, it generates 9 different numbers which correspond to 9 sheet impedance values in one supercell.

In the below [figure](https://github.com/Javadio/Discrete_sheet_impedance/blob/main/zvalues.jpg), you can see the generated 9 values for ${Z_{1...K}}$.
<h2>Figure: Z Values Visualization</h2>

<p align="center">
  <img src="https://github.com/Javadio/Discrete_sheet_impedance/blob/main/zvalues.jpg" alt="Z Values" width="900">
</p>

<p align="center"><b>Figure 1:</b> Visualization of Z values from the Discrete Sheet Impedance Model.</p>
---------------------------------------------------------------------------------------

Then by using check8GHz.mat code and by inserting the generated Z values as x matrix, you can calculate the amplitude of each scattered harmonics both propagating and evanescent waves. In the below [figure](https://github.com/Javadio/Discrete_sheet_impedance/blob/main/fig2.jpg), (a) and (c) you can see the generated harmonic amplitudes for 20 harmonics. Also (b) and (d) shows the results generated in CST.   

<h2>Figure: Harmonic amplitude Visualization</h2>

<p align="center">
  <img src="https://github.com/Javadio/Discrete_sheet_impedance/blob/main/fig2.jpg" alt="Z Values" width="900">
</p>

<p align="center"><b>Figure 2:</b> Visualization of harmonic amplitude based on the generated Z values.</p>

# Files
-------------------------------------------------------------------------------------
optim8GHz.m contains information on the MultiStart optimization framework to generate required ${Z_{1...K}}$. check8GHz.m file contains information on the calculation of the harmonic amplitude scattered from the periodic metasurface.  

# References
-------------------------------------------------------------------------------------
[A] [Shabanpour, J., Simovski, C., & Geraci, G. (2024). Physically Consistent RIS: From Reradiation Mode Optimization to Practical Realization. arXiv preprint arXiv:2409.17738](<https://arxiv.org/abs/2409.17738>).

[B] Shabanpour, Javad, X. Wang, S. Kosulnikov, and C. Simovski, “Advancing RIS beamforming efficiency: Moving beyond diagonal matrix techniques”, IEEE Transactions on Vehicular Technology, 2024, Under review.
