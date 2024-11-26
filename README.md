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
In the file optim8GHz.mat, you can define your frequency, number of the unit cells in the supecell, desired deflection frequency, substrate thickness and all the required information. By defining the cost function as $F({Z_{1...K}}) = \left| {{A_{\rm cal}} - {A_{\rm goal}}} \right|$ and employing the MultiStart and Fmincon optimization algorithms, Matlab searches the proper matrix ${Z_{1...K}}$ that minimizes the cost function i.e. maximal reflected power  in the prescribed direction w.r.t the incident power. 
