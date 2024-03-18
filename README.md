# Code-BD-RIS-SPAWC2024
This is a code package related to the paper "MIMO capacity maximization with beyond-diagonal RIS"  I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutierrez, submitted to SPAWC 2024.

# Content of Code Package
The code is implemented in Matlab, you might need cvx (with all subfolders) to run the code.
The files Script_SPAWC24vsBDRISposition, Script_SPAWC24vsM, ScriptSPAWC243vsPower and ScriptSPAWC24_Convergence are scripts that reproduce some of the plots in the paper. The programs compare the Max-capacity BDRIS algorithm proposed in the paper with a single-connected RIS, optimized using the algorithm in S. Zhang, R. Zhang, "Capacity characterization for IRS aided MIMO Communications", JSAC, vol. 38, no. 8, pp. 1823-1838, 2020.

The main function to optimize the BD-RIS is MaxCap_BDRIS_passive.
The main function for optimizing the diagonal RIS is MaxCap_RIS_passive.

# License and Referencing
This code package is licensed under the GNUv3 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

# Acknowledments
This work is supported by the European Commission’s Horizon Europe, Smart Networks and Services Joint Undertaking, research and innovation program under grant agreement 101139282, 6G-SENSES project. The work of I. Santamaria was also partly supported under grant PID2022-137099NB-C43 (MADDIE) funded by MCIN/AEI /10.13039/501100011033.
