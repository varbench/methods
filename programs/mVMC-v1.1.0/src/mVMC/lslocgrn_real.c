/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

his program is developed based on the mVMC-mini program
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
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * local Q and Green Functions with Lanczos step
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "lslocgrn_real.h"
#ifndef _SRC_LSLOCGRN_REAL
#define _SRC_LSLOCGRN_REAL
#include "calham_real.h"
#include "locgrn_real.c"
#include "pfupdate_two_real.h"
#include <complex.h>
#include "global.h"
#include "locgrn_real.h"
#include "workspace.c"
#include "calham_real.c"
#include "pfupdate_real.c"
#include "pfupdate_two_real.c"
#include "qp_real.c"
#include "projection.h"

/* Calculate <psi|QQ|x>/<psi|x> */
void LSLocalQ_real(const double h1, const double ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, double *_LSLQ_real)
{
  double e0,h2;

  e0 = CalculateHamiltonian0_real(eleNum); /* V */

  h2 = h1*e0; /* HV = (V+K+W)V */
  h2 += calculateHK_real(h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  h2 += calculateHW_real(h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);

  /* calculate local Q (IQ) */
  _LSLQ_real[0] = 1.0; /* I */
  _LSLQ_real[1] = h1;  /* H = V+K+W */

  /* calculate local Q (KQ) */
  _LSLQ_real[2] = h1;  /* H */
  _LSLQ_real[3] = h2;  /* H*H */

  return;
}

///
/// \param h1
/// \param ip
/// \param eleIdx
/// \param eleCfg
/// \param eleNum
/// \param eleProjCnt
/// \return val
/// \version 1.0
double calculateHK_real(const double h1, const double ip, int *eleIdx, int *eleCfg,
                           int *eleNum, int *eleProjCnt) {
  int idx,ri,rj,s;
  double val=0.0;

  for(idx=0;idx<NTransfer;idx++) {
    ri = Transfer[idx][0];
    rj = Transfer[idx][2];
    s  = Transfer[idx][3];

    val -= creal(ParaTransfer[idx]) * calHCA_real(ri,rj,s,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    /* Caution: negative sign */
  }

  return val;
}

/// calculate <psi| H C_is A_js |x>/<psi|x>
/// \param ri
/// \param rj
/// \param s
/// \param h1
/// \param ip
/// \param eleIdx
/// \param eleCfg
/// \param eleNum
/// \param eleProjCnt
/// \return val
/// \version 1.0
double calHCA_real(const int ri, const int rj, const int s,
                      const double h1, const double ip, int *eleIdx, int *eleCfg,
                      int *eleNum, int *eleProjCnt) {
  int rsi=ri+s*Nsite;
  int rsj=rj+s*Nsite;
  double val;
  double g;

  /* check */
  if(rsi==rsj) {
    if(eleNum[rsi]==1) return h1;
    else return 0.0;
  } else {
    if(eleNum[rsj]==0) return 0.0;
    if(eleNum[rsi]==1) return 0.0;
  }

  g = checkGF1_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum);
  if(fabs(g)>1.0e-12) {
    val = calHCA1_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  } else {
    val = calHCA2_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  }

  return val;
}

///
/// \param ri
/// \param rj
/// \param s
/// \param ip
/// \param eleIdx
/// \param eleCfg
/// \param eleNum
/// \return z/ip
/// \version 1.0
double checkGF1_real(const int ri, const int rj, const int s, const double ip,
                int *eleIdx, const int *eleCfg, int *eleNum) {
  double z;
  int mj,msj,rsi,rsj;
  double pfMNew[NQPFull];

  mj = eleCfg[rj+s*Nsite];
  msj = mj + s*Ne;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;

  /* hopping */
  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;

  /* calculate Pfaffian */
  CalculateNewPfM_real(mj, s, pfMNew, eleIdx, 0, NQPFull);
  z = CalculateIP_real(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return z/ip;
}

/// calculate <psi| H C_is A_js |x>/<psi|x> = <psi|x'>/<psi|x> * <psi|H|x'>/<psi|x'>
/// \param ri
/// \param rj
/// \param s
/// \param ip
/// \param eleIdx
/// \param eleCfg
/// \param eleNum
/// \param eleProjCnt
/// \return
/// \version 1.0
double calHCA1_real(const int ri, const int rj, const int s,
               const double ip, int *eleIdx, int *eleCfg,
               int *eleNum, int *eleProjCnt) {
  double *oldInvM; /* [NQPFull*Nsize*Nsize;] */
  double *oldPfM;  /* [NQPFull] */
  int *projCntNew;

  int rsi=ri+s*Nsite;
  int rsj=rj+s*Nsite;
  int mj;
  double ipNew,z,e;

  RequestWorkSpaceInt(NProj);
  RequestWorkSpaceDouble(NQPFull*(Nsize*Nsize+1));

  projCntNew = GetWorkSpaceInt(NProj);
  oldInvM = GetWorkSpaceDouble(NQPFull*Nsize*Nsize);
  oldPfM  = GetWorkSpaceDouble(NQPFull);

  /* copy InvM and PfM */
  copyMAll_real(InvM_real,PfM_real,oldInvM,oldPfM);

  /* The mj-th electron with spin s hops to site ri */
  mj = eleCfg[rsj];
  eleIdx[mj+s*Ne] = ri;
  eleCfg[rsj] = -1;
  eleCfg[rsi] = mj;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;

  UpdateProjCnt(rj, ri, s, projCntNew, eleProjCnt, eleNum);
  z = ProjRatio(projCntNew,eleProjCnt);

  UpdateMAll_real(mj,s,eleIdx,0,NQPFull);
  ipNew = CalculateIP_real(PfM_real,0,NQPFull,MPI_COMM_SELF);

  e = CalculateHamiltonian_real(ipNew,eleIdx,eleCfg,eleNum,projCntNew);

  /* revert hopping */
  eleIdx[mj+s*Ne] = rj;
  eleCfg[rsj] = mj;
  eleCfg[rsi] = -1;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  /* restore InvM and PfM */
  copyMAll_real(oldInvM,oldPfM,InvM_real,PfM_real);

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceDouble();
  return e*z*ipNew/ip;
}

/* calculate <psi| H C_is A_js |x>/<psi|x> for <psi|CA|x>/<psi|x>=0 */
/* Assuming ri!=rj, eleNum[rsi]=1, eleNum[rsj]=0 */
double calHCA2_real(const int ri, const int rj, const int s,
                       const double ip, int *eleIdx, int *eleCfg,
                       int *eleNum, int *eleProjCnt) {
  const int nsize=Nsize;
  const int nsite2=Nsite2;

  int idx;
  int rk,rl,sk;
  //    double complex val=0.0;
  double val=0.0;
  int rsi = ri+s*Nsite;
  int rsj = rj+s*Nsite;

  double g;

  double *buffer;
  int *bufferInt;

  int *myEleIdx, *myEleNum, *myBufferInt, *myRsi, *myRsj;
  double *myBuffer;
  //    double complex myValue=0;
  //    double complex v=0.0;
  double myValue=0;
  double v=0.0;


  RequestWorkSpaceInt(NProj);      /* for GreenFunc1 */
  RequestWorkSpaceDouble(NQPFull); /* for GreenFunc1 */
  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj+6);
  RequestWorkSpaceThreadDouble(NQPFull+3*Nsize);

  bufferInt = GetWorkSpaceInt(NProj);
  buffer = GetWorkSpaceDouble(NQPFull);

  /* H0 term */
  /* <psi|H0 CA|x>/<psi|x> = H0(x') <psi|CA|x>/<psi|x> */
  g = GreenFunc1_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum,eleProjCnt,bufferInt,buffer);

  /* hopping */
  eleNum[rsi] = 1;
  eleNum[rsj] = 0;

  val = g*CalculateHamiltonian0_real(eleNum);

  /* revert hopping */
  eleNum[rsi] = 0;
  eleNum[rsj] = 1;

  /* end of H0 term */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myBufferInt,myBuffer,myValue,myRsi,myRsj)  \
  reduction(+:v)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myBufferInt = GetWorkSpaceThreadInt(NProj);
    myRsi = GetWorkSpaceThreadInt(3);
    myRsj = GetWorkSpaceThreadInt(3);
    myBuffer = GetWorkSpaceThreadDouble(NQPFull+3*Nsize);

#pragma loop noalias
    for(idx=0;idx<nsize;idx++) myEleIdx[idx] = eleIdx[idx];
#pragma loop noalias
    for(idx=0;idx<nsite2;idx++) myEleNum[idx] = eleNum[idx];

    myValue = 0.0;

    /* Transfer */
#pragma omp for private(idx,rk,rl,sk) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      rk = Transfer[idx][0];
      rl = Transfer[idx][2];
      sk = Transfer[idx][3];

      myValue -= creal(ParaTransfer[idx])
                       * GreenFunc2_real(rk,rl,ri,rj,sk,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBufferInt,myBuffer);
      /* Caution: negative sign */
    }

    /* Pair Hopping */
#pragma omp for private(idx,rk,rl) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      rk = PairHopping[idx][0];
      rl = PairHopping[idx][1];
      myRsi[0] = rk; /* s=0 */
      myRsj[0] = rl; /* s=0 */
      myRsi[1] = rk+Nsite; /* s=1 */
      myRsj[1] = rl+Nsite; /* s=1 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;

      myValue += ParaPairHopping[idx]
        * GreenFuncN_real(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);
    }

    /* Exchange Coupling */
#pragma omp for private(idx,rk,rl) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      rk = ExchangeCoupling[idx][0];
      rl = ExchangeCoupling[idx][1];
      myRsi[0] = rk; /* s=0 */
      myRsj[0] = rl; /* s=0 */
      myRsi[1] = rl+Nsite; /* s=1 */
      myRsj[1] = rk+Nsite; /* s=1 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myValue += ParaExchangeCoupling[idx]
        * GreenFuncN_real(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);

      myRsi[0] = rk+Nsite; /* s=1 */
      myRsj[0] = rl+Nsite; /* s=1 */
      myRsi[1] = rl; /* s=0 */
      myRsj[1] = rk; /* s=0 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myValue += ParaExchangeCoupling[idx]
        * GreenFuncN_real(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);
    }

    /* Inter All */
#pragma omp for private(idx) schedule(dynamic) nowait
    for(idx=0;idx<NInterAll;idx++) {
      myRsi[0] = InterAll[idx][0] + InterAll[idx][3]*Nsite;
      myRsj[0] = InterAll[idx][2] + InterAll[idx][3]*Nsite;
      myRsi[1] = InterAll[idx][4] + InterAll[idx][7]*Nsite;
      myRsj[1] = InterAll[idx][6] + InterAll[idx][7]*Nsite;
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myValue += creal(ParaInterAll[idx])
        * GreenFuncN_real(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);
    }

    v += myValue;
  }
  val += v;

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceDouble();
  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();

  return val;
}


/* copy invM and pfM */
void copyMAll_real(double *invM_from, double *pfM_from, double *invM_to, double *pfM_to) {
  int i,n;

  n = NQPFull*Nsize*Nsize;
#pragma loop noalias
  for(i=0;i<n;i++) invM_to[i] = invM_from[i];

  n = NQPFull;
#pragma loop noalias
  for(i=0;i<n;i++) pfM_to[i] = pfM_from[i];

  return;
}

double calculateHW_real(const double h1, const double ip, int *eleIdx, int *eleCfg,
                           int *eleNum, int *eleProjCnt) {
  int idx,ri,rj,s,rk,rl,t;
  double val=0.0,tmp;

  /* Pair Hopping */
  for(idx=0;idx<NPairHopping;idx++) {
    ri = PairHopping[idx][0];
    rj = PairHopping[idx][1];

    val += ParaPairHopping[idx]
      * calHCACA_real(ri,rj,ri,rj,0,1,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  }

  /* Exchange Coupling */
  for(idx=0;idx<NExchangeCoupling;idx++) {
    ri = ExchangeCoupling[idx][0];
    rj = ExchangeCoupling[idx][1];

    tmp =  calHCACA_real(ri,rj,rj,ri,0,1,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    tmp += calHCACA_real(ri,rj,rj,ri,1,0,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    val += ParaExchangeCoupling[idx] * tmp;
  }

  /* Inter All */
  for(idx=0;idx<NInterAll;idx++) {
    ri = InterAll[idx][0];
    rj = InterAll[idx][2];
    s  = InterAll[idx][3];
    rk = InterAll[idx][4];
    rl = InterAll[idx][6];
    t  = InterAll[idx][7];

    val += creal(ParaInterAll[idx])
      * calHCACA_real(ri,rj,rk,rl,s,t,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  }

  return val;
}


double calHCACA_real(const int ri, const int rj, const int rk, const int rl,
                        const int si,const int sk,
                        const double h1, const double ip, int *eleIdx, int *eleCfg,
                        int *eleNum, int *eleProjCnt) {
  int rsi=ri+si*Nsite;
  int rsj=rj+si*Nsite;
  int rsk=rk+sk*Nsite;
  int rsl=rl+sk*Nsite;

  double val;
  double g;

  /* check */
  if(rsk==rsl) {
    if(eleNum[rsk]==1) {
      return calHCA_real(ri,rj,si,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    } else return 0;
  } else if(rsj==rsk) {
    if(eleNum[rsj]==1) return 0;
    else {
      return calHCA_real(ri,rl,si,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    }
  } else if(rsj==rsl) {
    return 0;
  } else if(rsi==rsj) {
    if(eleNum[rsi]==1) {
      return calHCA_real(rk,rl,sk,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    } else return 0;
  } else if(rsi==rsk) {
    return 0;
  } else if(rsi==rsl) {
    if(eleNum[rsi]==1) {
      return -calHCA_real(rk,rj,sk,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    } else return 0;
  } else {
    if(eleNum[rsl]==0) return 0.0;
    if(eleNum[rsk]==1) return 0.0;
    if(eleNum[rsj]==0) return 0.0;
    if(eleNum[rsi]==1) return 0.0;
  }

  g = checkGF2_real(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum);
  if(fabs(g)>1.0e-12) {
    val = calHCACA1_real(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  } else {
    val = calHCACA2_real(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  }

  return val;
}

double checkGF2_real(const int ri, const int rj, const int rk, const int rl,
                const int s, const int t, const double ip,
                int *eleIdx, const int *eleCfg, int *eleNum) {
  double z;
  int mj,msj,ml,mtl;
  int rsi,rsj,rtk,rtl;
  double *pfMNew;
  double *buffer;

  RequestWorkSpaceDouble(NQPFull+2*Nsize);
  pfMNew = GetWorkSpaceDouble(NQPFull);
  buffer = GetWorkSpaceDouble(2*Nsize);

  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;
  rtk = rk + t*Nsite;
  rtl = rl + t*Nsite;

  mj = eleCfg[rj+s*Nsite];
  ml = eleCfg[rl+t*Nsite];
  msj = mj + s*Ne;
  mtl = ml + t*Ne;

  /* hopping */
  eleIdx[mtl] = rk;
  eleNum[rtl] = 0;
  eleNum[rtk] = 1;

  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;

  /* calculate Pfaffian */
  CalculateNewPfMTwo_real(ml, t, mj, s, pfMNew, eleIdx, 0, NQPFull, buffer);
  z = CalculateIP_real(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[mtl] = rl;
  eleNum[rtl] = 1;
  eleNum[rtk] = 0;
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  ReleaseWorkSpaceDouble();
  return z/ip;
}

double calHCACA1_real(const int ri, const int rj, const int rk, const int rl,
                 const int si,const int sk,
                 const double ip, int *eleIdx, int *eleCfg,
                 int *eleNum, int *eleProjCnt) {
  double *oldInvM; /* [NQPFull*Nsize*Nsize;] */
  double *oldPfM;  /* [NQPFull] */
  int *projCntNew;

  int rsi=ri+si*Nsite;
  int rsj=rj+si*Nsite;
  int rsk=rk+sk*Nsite;
  int rsl=rl+sk*Nsite;
  int mj,ml;
  double ipNew,z,e;

  RequestWorkSpaceInt(NProj);
  RequestWorkSpaceDouble(NQPFull*(Nsize*Nsize+1));

  projCntNew = GetWorkSpaceInt(NProj);
  oldInvM = GetWorkSpaceDouble(NQPFull*Nsize*Nsize);
  oldPfM  = GetWorkSpaceDouble(NQPFull);

  /* copy InvM and PfM */
  copyMAll_real(InvM_real,PfM_real,oldInvM,oldPfM);

  /* The ml-th electron with spin sk hops from rl to rk */
  ml = eleCfg[rsl];
  eleIdx[ml+sk*Ne] = rk;
  eleCfg[rsl] = -1;
  eleCfg[rsk] = ml;
  eleNum[rsl] = 0;
  eleNum[rsk] = 1;
  UpdateProjCnt(rl, rk, sk, projCntNew, eleProjCnt, eleNum);

  /* The mj-th electron with spin si hops from rj to ri */
  mj = eleCfg[rsj];
  eleIdx[mj+si*Ne] = ri;
  eleCfg[rsj] = -1;
  eleCfg[rsi] = mj;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, si, projCntNew, projCntNew, eleNum);

  z = ProjRatio(projCntNew,eleProjCnt);

  UpdateMAllTwo_real(ml, sk, mj, si, rl, rj, eleIdx, 0, NQPFull);
  ipNew = CalculateIP_real(PfM_real,0,NQPFull,MPI_COMM_SELF);

  e = CalculateHamiltonian_real(ipNew,eleIdx,eleCfg,eleNum,projCntNew);

  /* revert hopping */
  eleIdx[mj+si*Ne] = rj;
  eleCfg[rsj] = mj;
  eleCfg[rsi] = -1;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  eleIdx[ml+sk*Ne] = rl;
  eleCfg[rsl] = ml;
  eleCfg[rsk] = -1;
  eleNum[rsl] = 1;
  eleNum[rsk] = 0;

  /* restore InvM and PfM */
  copyMAll_real(oldInvM,oldPfM,InvM_real,PfM_real);

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceDouble();
  return e*z*ipNew/ip;
}

/* calculate <psi| H C_is A_js C_kt A_lt|x>/<psi|x> for <psi|CACA|x>/<psi|x>=0 */
/* Assuming ri,rj,rk,rl are different, eleNum[rsi]=1, eleNum[rsj]=0, eleNum[rsk]=1, eleNum[rsl]=0  */
double calHCACA2_real(const int ri, const int rj, const int rk, const int rl,
                         const int si,const int sk,
                         const double ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt) {
  const int nsize=Nsize;
  const int nsite2=Nsite2;

  const int rsi = ri+si*Nsite;
  const int rsj = rj+si*Nsite;
  const int rsk = rk+sk*Nsite;
  const int rsl = rl+sk*Nsite;

  int idx,r0,r1;
  double val=0.0;

  double g;

  double *buffer;
  int *bufferInt;

  int *myEleIdx, *myEleNum, *myBufferInt, *myRsi, *myRsj;
  double *myBuffer;
  double myValue=0.0;
  double v=0.0;

  RequestWorkSpaceInt(NProj);      /* for GreenFunc2 */
  RequestWorkSpaceDouble(NQPFull+2*Nsize); /* for GreenFunc2 */
  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj+8);
  RequestWorkSpaceThreadDouble(NQPFull+3*Nsize);

  bufferInt = GetWorkSpaceInt(NProj);
  buffer = GetWorkSpaceDouble(NQPFull+2*Nsize);

  /* H0 term */
  /* <psi|H0 CACA|x>/<psi|x> = H0(x') <psi|CACA|x>/<psi|x> */
  g = GreenFunc2_real(ri,rj,rk,rl,si,sk,ip,
      eleIdx,eleCfg,eleNum,eleProjCnt,bufferInt,buffer);

  /* hopping */
  eleNum[rsi] = 1;
  eleNum[rsj] = 0;
  eleNum[rsk] = 1;
  eleNum[rsl] = 0;

  val = g*CalculateHamiltonian0_real(eleNum);

  /* revert hopping */
  eleNum[rsi] = 0;
  eleNum[rsj] = 1;
  eleNum[rsk] = 0;
  eleNum[rsl] = 1;

  /* end of H0 term */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myBufferInt,myBuffer,myValue,myRsi,myRsj)  \
  reduction(+:v)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myBufferInt = GetWorkSpaceThreadInt(NProj);
    myRsi = GetWorkSpaceThreadInt(4);
    myRsj = GetWorkSpaceThreadInt(4);
    myBuffer = GetWorkSpaceThreadDouble(NQPFull+4*Nsize);

#pragma loop noalias
    for(idx=0;idx<nsize;idx++) myEleIdx[idx] = eleIdx[idx];
#pragma loop noalias
    for(idx=0;idx<nsite2;idx++) myEleNum[idx] = eleNum[idx];

    myValue = 0.0;

    /* Transfer */
#pragma omp for private(idx) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      myRsi[0] = Transfer[idx][0]+Transfer[idx][1]*Nsite;
      myRsj[0] = Transfer[idx][2]+Transfer[idx][3]*Nsite;
      myRsi[1] = rsi;
      myRsj[1] = rsj;
      myRsi[2] = rsk;
      myRsj[2] = rsl;

      myValue -= creal(ParaTransfer[idx])
        * GreenFuncN_real(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);
      /* Caution: negative sign */
    }

    /* Pair Hopping */
#pragma omp for private(idx,r0,r1) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      r0 = PairHopping[idx][0];
      r1 = PairHopping[idx][1];
      myRsi[0] = r0; /* s=0 */
      myRsj[0] = r1; /* s=0 */
      myRsi[1] = r0+Nsite; /* s=1 */
      myRsj[1] = r1+Nsite; /* s=1 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myRsi[3] = rsk;
      myRsj[3] = rsl;

      myValue += ParaPairHopping[idx]
        * GreenFuncN_real(4,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);
    }

    /* Exchange Coupling */
#pragma omp for private(idx,r0,r1) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      r0 = ExchangeCoupling[idx][0];
      r1 = ExchangeCoupling[idx][1];
      myRsi[0] = r0; /* s=0 */
      myRsj[0] = r1; /* s=0 */
      myRsi[1] = r1+Nsite; /* s=1 */
      myRsj[1] = r0+Nsite; /* s=1 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myRsi[3] = rsk;
      myRsj[3] = rsl;
      myValue += ParaExchangeCoupling[idx]
        * GreenFuncN_real(4,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);

      myRsi[0] = r0+Nsite; /* s=1 */
      myRsj[0] = r1+Nsite; /* s=1 */
      myRsi[1] = r1; /* s=0 */
      myRsj[1] = r0; /* s=0 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myRsi[3] = rsk;
      myRsj[3] = rsl;
      myValue += ParaExchangeCoupling[idx]
        * GreenFuncN_real(4,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);
    }

    /* Inter All */
#pragma omp for private(idx) schedule(dynamic) nowait
    for(idx=0;idx<NInterAll;idx++) {
      myRsi[0] = InterAll[idx][0] + InterAll[idx][3]*Nsite;
      myRsj[0] = InterAll[idx][2] + InterAll[idx][3]*Nsite;
      myRsi[1] = InterAll[idx][4] + InterAll[idx][7]*Nsite;
      myRsj[1] = InterAll[idx][6] + InterAll[idx][7]*Nsite;
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myRsi[3] = rsk;
      myRsj[3] = rsl;
      myValue += creal(ParaInterAll[idx])
        * GreenFuncN_real(4,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);
    }

    v += myValue;
  }
  val += v;

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceDouble();
  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();

  return val;
}


/* Calculate <psi|QCisAjs|x>/<psi|x> */
void LSLocalCisAjs_real(const double h1, const double ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt) {
  const int nCisAjs=NCisAjs;
  double *lsLCisAjs_real = LSLCisAjs_real;
  double complex*localCisAjs = LocalCisAjs;
  int ri,rj,s;
  int idx;

  /* copy local ICisAjs */
#pragma loop noalias
  for(idx=0;idx<nCisAjs;idx++){
    lsLCisAjs_real[idx] = creal(localCisAjs[idx]);
  }

  for(idx=0;idx<nCisAjs;idx++){
    ri = CisAjsIdx[idx][0];
    rj = CisAjsIdx[idx][2];
    s  = CisAjsIdx[idx][3];

    /* calculate local HCisAjs */
    LSLCisAjs_real[idx+nCisAjs] = calHCA_real(ri,rj,s,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  }
  return;
}
#endif
