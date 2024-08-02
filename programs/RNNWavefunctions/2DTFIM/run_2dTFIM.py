from Training2DRNN_2DTFIM import run_2DTFIM

#numsteps = number of training iterations
#systemsize_x = the size of the x-dimension of the square lattice
#systemsize_x = the size of the y-dimension of the square lattice
#Bx = transverse magnetic field
#numsamples = number of samples used for training
#num_units = number of memory units of the hidden state of the RNN
#num_layers is not supported yet, stay tuned!
RNNEnergy, varRNNEnergy, meanE_final, varEnergy  = run_2DTFIM(numsteps = 150000, systemsize_x = 12, systemsize_y = 12, Bx = 3, num_units = 500, numsamples = 500, learningrate = 5e-3, seed = 111)

#RNNEnergy is a numpy array of the variational energy of the pRNN wavefunction
#varRNNEnergy is a numpy array of the variance of the variational energy of the pRNN wavefunction
#meanE_final is the final approximation of the ground state energy after training
#varEnergy is the final estimation of the energy variance after training
