afqmc_param_*dt.txt contain the input parameters of AFQMC. * in the file name refers to the imaginary time steps dt used in calculations. The results are used to extrapolate to dt=0. The variable names are listed in file afqmc_param_names.txt

phi.dat is the initial wavefunction. The first line is the logrithm of its weight, where the first column is the real part and the second column is the imaginary part. In this case the weight is 1 because we only have a single Slater determinant in the initial wavefunction. The second line means that the dimension of the data is two-dimensional. The first column in the third line describes the number of elements in each electron orbital (number of sites * spin).  The second column in the third line is the number of electrons. The wavefunction starts at the fourth line and is written as one orbital after another.

Please request permission from the original author for the code at https://github.com/hshi/AFQMCLAB
