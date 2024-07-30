/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.
*/
#include "initial.h"
#include "../../sfmt/SFMT.h"

void initial(struct BindStruct *X) {

  int int_i, int_j;
  int spin_0, spin_1;
  int site_0, site_1;
  int t_site_0, t_site_1;
  int Ns;
  double complex tmp;

  Ns = X->Def.Nsite;

  if (X->Def.NInitial == 0) {
    for (int_i = 0; int_i < 2 * X->Def.Nsite; int_i++) {
      for (int_j = 0; int_j < 2 * X->Def.Nsite; int_j++) {
        //X->Large.G[int_i][int_j]   = 0.0;
        X->Large.G[int_i][int_j] = 0.01 * (genrand_real2() - 0.5); /* uniform distribution [0,1) */
        //X->Large.G[int_i][int_j]  += 1*I*genrand_real2(); /* uniform distribution [0,1) */
        //X->Large.G[int_i][int_j]  /=sqrt(2.0);
      }
    }
  }

  if(X->Def.print==1){
    printf("#[s]output initial Green functions\n");
  }
  //Case: Initial green's functions are defined.
  for (int_i = 0; int_i < X->Def.NInitial; int_i++) {
    site_0 = X->Def.Initial[int_i][0];
    spin_0 = X->Def.Initial[int_i][1];
    site_1 = X->Def.Initial[int_i][2];
    spin_1 = X->Def.Initial[int_i][3];
    tmp = X->Def.ParaInitial[int_i];

    if(X->Def.print==1){
      printf("%d %lf %lf \n", int_i, creal(tmp), cimag(tmp));
    }
    t_site_0 = site_0 + spin_0 * Ns;
    t_site_1 = site_1 + spin_1 * Ns;
    X->Large.G[t_site_0][t_site_1] = tmp;
  }
  if(X->Def.print==1){
    printf("#[e]output initial Green functions\n");
  }
}
