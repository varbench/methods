This is a script for calculating the ground-state energy of the simple Hubbard model using mVMC.
For simplicity, we focus on the case of Ns=8x8, Ne=50 or 64, U=4 or 8 (the data in https://github.com/varbench/varbench/blob/main/Hubbard/supplements/Imada_group_overview.md ).
Please run the python code make_all.py to generate definition files.
The code requires numpy and json libraries.
Input parameters are written in dat_input_makedef.txt.

Here is a tip for optimization:
1. Optimize f_ij of pair wave function only using the unrestricted Hartree-Fock solution.
2. Once you get the optimized f_ij, turn on RBM parameters and optimize them.
3. Once you get the optimized f_ij and RBM parameters, turn on the spin and momentum projections and optimize f_ij and RBM again.
4. When you optimize all the parameters, apply the Lanczos method to calculate the ground-state energy.
