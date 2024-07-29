c---------------------------------------------------------------------
c The true Deltas are twice the ones used here!
c---------------------------------------------------------------------
      program prepv
      implicit none

c *** Lattice
      INTEGER(4) nh,nbase,ndim
      INTEGER(4) nindm,nindms,maxmulti,maxmultis
      REAL(8) pi,argo
      INTEGER(4), dimension(:), allocatable :: imulti,imultis
      INTEGER(4), dimension(:,:,:), allocatable :: ivicn,ivicns
      REAL(8), dimension(:), allocatable :: x
      REAL(8), dimension(:), allocatable :: qxv
      COMPLEX(8), dimension(:,:,:), allocatable :: phaserKIN,phaserBCS
      COMPLEX(8), dimension(:,:,:), allocatable :: phaseiKIN,phaseiBCS
c *** State
      INTEGER(4) nelup,neldo
      INTEGER(4) iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf
      REAL(8) muf,Qx
      REAL(8), dimension(:), allocatable :: phi
      REAL(8), dimension(:), allocatable :: tvar
      REAL(8), dimension(:), allocatable :: dBCS_r,dek_r,donsite_r,dmu
      REAL(8), dimension(:), allocatable :: dBCS_i,dek_i,donsite_i
      REAL(8), dimension(:), allocatable :: dAF
      REAL(8), dimension(:), allocatable :: vjspin,vjdens
      REAL(8), dimension(:), allocatable :: eig
      REAL(8), dimension(:,:), allocatable :: v,vsz
      COMPLEX(8), dimension(:,:), allocatable :: z
      COMPLEX(8), dimension(:), allocatable :: alphai
c *** Hamiltonian
      INTEGER(4) irange,ioptPH
      REAL(8) uvv,jperp
      REAL(8), dimension(:), allocatable :: tham,jham
c *** Other
      INTEGER(4) i,j,k,icheck,icount
      REAL(8) rnup,rndo,dnrm2

      namelist /rangeham/ irange,ioptPH
      namelist /couplings/ jperp,uvv
      namelist /nelectrons/ nelup,neldo
      namelist /wavefunction/ iesspin,iesdens,iesBCS,ieskin,iesmu
     >,iesonsite,iesaf,muf

c*****************************************************************************
c BEGIN READ FROM FILES
c*****************************************************************************
      open(unit=10,file='fort.10',form='unformatted',status='unknown') 

c     ndim is the number of sites in the Bravais lattice
c     nbase is the number of sites in the basis
c     nh is the TOTAL number of sites
c     nindm is the number of GEOMETRICAL distances (excluding 0)
c     nindms is the number of distances with sign, calculated in LATTICE.f
c     irange is the number of couplings in the Hamiltonian
      read(9) ndim,nbase,nindm,nindms,maxmulti,maxmultis

      write(6,*) 'Number of sites in the Bravais lattice=',ndim
      write(6,*) 'Number of sites in the basis          =',nbase
      write(6,*) 'Total number of sites                 =',ndim*nbase

      nh=ndim*nbase

      ALLOCATE(imulti(nindm))
      ALLOCATE(imultis(nindms))
      ALLOCATE(ivicn(nh,maxmulti,nindm))
      ALLOCATE(ivicns(nh,maxmultis,nindms))
      ALLOCATE(phaserKIN(nh,maxmultis,nindms))
      ALLOCATE(phaseiKIN(nh,maxmultis,nindms))
      ALLOCATE(phaserBCS(nh,maxmultis,nindms))
      ALLOCATE(phaseiBCS(nh,maxmultis,nindms))

      read(9) imulti,imultis,ivicn,ivicns
      read(9) phaserKIN,phaserBCS,phaseiKIN,phaseiBCS

      ALLOCATE(x(ndim))
      ALLOCATE(qxv(ndim))
      ALLOCATE(phi(nbase))

      do i=1,ndim
       read(27,*) qxv(i)
       read(28,*) x(i)
      enddo

      ioptPH=0           ! ioptPH.ge.0 PH ioptPH.lt.0 NO PH
      read(5,rangeham)

      ALLOCATE(tham(irange))
      ALLOCATE(jham(irange))
      ALLOCATE(tvar(irange))

      nelup=0
      neldo=0
      iesspin=0
      iesdens=0
      iesBCS=0
      ieskin=0
      iesmu=0
      iesonsite=0
      iesaf=0
      muf=0.d0
      do i=1,irange
       tham(i)=0.d0
       jham(i)=0.d0
      enddo
      jperp=0.d0
      uvv=0.d0

      read(5,couplings)
      read(5,nelectrons)
      read(5,wavefunction)

      if(ioptPH.lt.0) then
       if(iesBCS.gt.0.or.iesonsite.gt.0) then
        write(6,*) 'No BCS terms with no Particle-hole!'
        stop
       endif
      elseif(ioptPH.ge.0) then
       if(iesaf.gt.0) then
        write(6,*) 'No AF terms with Particle-hole!'
        stop
       endif
      endif

      ALLOCATE(dBCS_r(iesBCS))
      ALLOCATE(dBCS_i(iesBCS))
      ALLOCATE(dek_r(ieskin))
      ALLOCATE(dek_i(ieskin))
      ALLOCATE(dmu(iesmu))
      ALLOCATE(donsite_r(iesonsite))
      ALLOCATE(donsite_i(iesonsite))
      ALLOCATE(dAF(iesaf))
      ALLOCATE(vjspin(iesspin))
      ALLOCATE(vjdens(iesdens))

      if(iesspin.gt.nindm) then
       write(6,*) 'iesspin too large',iesspin,nindm
       stop
      endif
      if(iesdens.gt.nindm) then
       write(6,*) 'iesdens too large',iesdens,nindm
       stop
      endif
      if(iesBCS.gt.nindms) then
       write(6,*) 'iesBCS too large',iesBCS,nindms
       stop
      endif
      if(ieskin.gt.nindms) then
       write(6,*) 'ieskin too large',ieskin,nindms
       stop
      endif
      if(iesmu.gt.1) then
       write(6,*) 'iesmu too large, it must be 0 or 1',iesmu
       stop
      endif
      if(iesonsite.gt.1) then
       write(6,*) 'iesonsite too large, it must be 0 or 1',iesonsite
       stop
      endif
      if(iesaf.gt.1) then
       write(6,*) 'iesaf too large, it must be 0 or 1',iesaf
       stop
      endif

      if(iesspin.ne.0) then 
       read(5,*) (vjspin(k),k=1,iesspin)   ! Jastrow parameters Spin 
      else
       read(5,*)
      endif
      if(iesdens.ne.0) then 
       read(5,*) (vjdens(k),k=1,iesdens)   ! Jastrow parameters Density
      else
       read(5,*)
      endif
      if(iesBCS.ne.0) then 
       read(5,*) (dBCS_r(k),k=1,iesBCS) ! BCS parameters
      else
       read(5,*)
      endif
      if(iesBCS.ne.0) then 
       read(5,*) (dBCS_i(k),k=1,iesBCS) ! BCS parameters
      else
       read(5,*)
      endif
      if(ieskin.ne.0) then
       read(5,*) (dek_r(k),k=1,ieskin) ! long-range hopping 
      else
       read(5,*)
      endif
      if(ieskin.ne.0) then
       read(5,*) (dek_i(k),k=1,ieskin) ! long-range hopping 
      else
       read(5,*)
      endif
      if(iesmu.ne.0) then
       read(5,*) (dmu(k),k=1,iesmu) ! chemical potential
      else
       read(5,*)
      endif
      if(iesonsite.ne.0) then
       read(5,*) (donsite_r(k),k=1,iesonsite) ! on-site pairing
      else
       read(5,*)
      endif
      if(iesonsite.ne.0) then
       read(5,*) (donsite_i(k),k=1,iesonsite) ! on-site pairing
      else
       read(5,*)
      endif
      if(iesaf.ne.0) then
       read(5,*) (dAF(k),k=1,iesaf) ! AF field
       read(5,*) Qx,(phi(k),k=1,nbase)
      else
       read(5,*)
       read(5,*)
      endif
      read(5,*) (tvar(k),k=1,irange)
      read(5,*) (tham(k),k=1,irange)
      read(5,*) (jham(k),k=1,irange)

      icheck=0
      do i=1,ndim
       if(dabs(Qx-qxv(i)).lt.1.d-7)
     >  icheck=icheck+1
      enddo
      if(icheck.ne.1) then
       write(6,*) 'Wrong momentum selected in geometry'
       stop
      endif

      ALLOCATE(alphai(nh))
      pi=dacos(-1.d0)
      icount=0
      do i=1,ndim
       do j=1,nbase
        icount=icount+1
        argo=2.d0*pi*(x(i)*Qx+phi(j))
        alphai(icount)=dcmplx(dcos(argo),dsin(argo))
       enddo
      enddo

      write(6,*)
      write(6,*) 'wavefunction parameters'
      write(6,*)

      write(6,*) (vjspin(k),k=1,iesspin)
      write(6,*) (vjdens(k),k=1,iesdens)
      write(6,*) (dBCS_r(k),k=1,iesBCS)
      write(6,*) (dBCS_i(k),k=1,iesBCS)
      write(6,*) (dek_r(k),k=1,ieskin)
      write(6,*) (dek_i(k),k=1,ieskin)
      write(6,*) (dmu(k),k=1,iesmu)
      write(6,*) (donsite_r(k),k=1,iesonsite)
      write(6,*) (donsite_i(k),k=1,iesonsite)
      write(6,*) (dAF(k),k=1,iesaf)
c*****************************************************************************
c END READ FROM FILES
c*****************************************************************************

      ALLOCATE(v(nh,nh))
      ALLOCATE(vsz(nh,nh))
      ALLOCATE(z(2*nh,2*nh))
      ALLOCATE(eig(2*nh))

      call ekdk(nh,irange,ivicn,ivicns,maxmulti,maxmultis
     >,imulti,imultis,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS
     >,iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf
     >,vjspin,vjdens,dBCS_r,dBCS_i,dek_r,dek_i,dmu,donsite_r
     >,donsite_i,dAF,alphai,vsz,v,z,muf,tvar,eig,ioptPH)

      write(6,*) '   i    Eigenvalue       Norm up         Norm down'
      do i=1,2*nh
       rnup=0.d0
       do j=1,nh
        rnup=rnup+z(j,i)*dconjg(z(j,i))
       enddo
       rnup=dsqrt(rnup)
       if(rnup.le.1d-10) rnup=0.d0
       rndo=0.d0
       do j=nh+1,2*nh
        rndo=rndo+z(j,i)*dconjg(z(j,i))
       enddo
       rndo=dsqrt(rndo)
       if(rndo.le.1d-10) rndo=0.d0
       write(6,999) i,eig(i),rnup,rndo
      enddo
999   format(i5,3x,3(f14.10,2x))

      rewind(10)
      write(10) ndim,nbase,nh,nindm,nindms,maxmulti,maxmultis
      write(10) ioptPH,irange,nelup,neldo
      write(10) imulti,imultis,ivicn,ivicns
      write(10) phaserKIN,phaserBCS,phaseiKIN,phaseiBCS
      write(10) iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf,muf
      write(10) z,eig,vjspin,vjdens,v,vsz,dBCS_r,dBCS_i,dek_r,dek_i,dmu
     >,donsite_r,donsite_i,dAF,alphai
      write(10) tham,jham,tvar,jperp,uvv

      close(10) 

      DEALLOCATE(imulti,imultis)
      DEALLOCATE(ivicn,ivicns)
      DEALLOCATE(phaserKIN,phaserBCS,phaseiKIN,phaseiBCS)
      DEALLOCATE(x,qxv,phi,alphai)
      DEALLOCATE(tham,jham,tvar)
      DEALLOCATE(dBCS_r,dBCS_i,dek_r,dek_i,dmu,donsite_r,donsite_i,dAF)
      DEALLOCATE(vjspin,vjdens)
      DEALLOCATE(z,eig,v,vsz)

      stop 
      end

      subroutine ekdk(nh,irange,ivicn,ivicns,maxmulti,maxmultis
     >,imulti,imultis,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS
     >,iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf
     >,vjspin,vjdens,dBCS_r,dBCS_i,dek_r,dek_i,dmu,donsite_r
     >,donsite_i,dAF,alphai,vsz,v,z,muf,tvar,eig,ioptPH)
      implicit none 
      integer*4 nh,iq,ip,info,i,k,j,jn,n1,maxmulti,maxmultis
     >,irange,ioptPH
      integer*4 iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf
      integer*4 ivicn(nh,maxmulti,*),ivicns(nh,maxmultis,*)
     >,imulti(*),imultis(*)
      real*8 vjdens(*),vjspin(*),dBCS_r(*),dBCS_i(*),dek_r(*),dek_i(*)
     >,dmu(*),donsite_r(*),donsite_i(*),dAF(*),tvar(*)
      real*8 vsz(nh,nh),v(nh,nh),eig(*),muf,dek1
      real*8 test1,test2
      complex(8) ai
      complex(8) alphai(nh)
      complex(8) z(2*nh,*)

      COMPLEX(8) phaserKIN(nh,maxmultis,*),phaserBCS(nh,maxmultis,*)
      COMPLEX(8) phaseiKIN(nh,maxmultis,*),phaseiBCS(nh,maxmultis,*)

      REAL(8), dimension(:), allocatable :: rpsip
      COMPLEX(8), dimension(:), allocatable :: psip
      COMPLEX(8), dimension(:,:), allocatable :: matrix

c Prepare the Jastrow terms

      do i=1,nh
       do j=1,nh
        vsz(i,j)=1.d0
        v(i,j)=1.d0
       enddo
      enddo
       
      do i=1,iesdens
       do k=1,imulti(i)
        do j=1,nh
         jn=abs(ivicn(j,k,i))
         if(jn.ne.0) then
          v(jn,j)=dexp(vjdens(i))
         endif
        enddo
       enddo
      enddo

      do i=1,iesspin
       do k=1,imulti(i)
        do j=1,nh
         jn=abs(ivicn(j,k,i))
         if(jn.ne.0) then 
          vsz(jn,j)=dexp(vjspin(i))
         endif
        enddo
       enddo
      enddo

c Prepare the Slater determinant

      ai=dcmplx(0.d0,1.d0)

      ALLOCATE(matrix(2*nh,2*nh))

      write(6,*)
      write(6,*) ' Chosen mu       =',muf
      if(iesmu.ge.1) then
       write(6,*) ' Renormalized mu =',muf+dmu(1)
      endif
      write(6,*)
         
      do ip=1,2*nh
       do iq=1,2*nh
        matrix(ip,iq)=dcmplx(0.d0,0.d0)
       enddo
      enddo

      do ip=1,nh
       matrix(ip,ip)=matrix(ip,ip)-muf
       if(ioptPH.ge.0) then
        matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)+muf
       elseif(ioptPH.lt.0) then
        matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)-muf
       endif
      enddo

      dek1=-1.d0
      do iq=1,irange
       do k=1,imulti(iq)
        do ip=1,nh
         jn=ivicn(ip,k,iq)
         if(jn.ne.0) then
          if(jn.gt.0) then 
           matrix(ip,jn)=matrix(ip,jn)+tvar(iq)*dek1
           if(ioptPH.ge.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)-tvar(iq)*dek1
           elseif(ioptPH.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)+tvar(iq)*dek1
           endif
          elseif(jn.lt.0) then
           jn=-jn
           matrix(ip,jn)=matrix(ip,jn)-tvar(iq)*dek1
           if(ioptPH.ge.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)+tvar(iq)*dek1
           elseif(ioptPH.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)-tvar(iq)*dek1
           endif
          endif
         endif
        enddo
       enddo
      enddo

      if(iesmu.ge.1) then
       do ip=1,nh
        matrix(ip,ip)=matrix(ip,ip)-dmu(1)
        if(ioptPH.ge.0) then
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)+dmu(1)
        elseif(ioptPH.lt.0) then
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)-dmu(1)
        endif
       enddo
      endif

      do iq=1,ieskin
       do k=1,imultis(iq)
        do ip=1,nh
         jn=ivicns(ip,k,iq)
         if(jn.ne.0) then
          if(jn.gt.0) then
           matrix(ip,jn)=matrix(ip,jn)
     >     +dek_r(iq)*phaserKIN(ip,k,iq)
     >     +ai*dek_i(iq)*phaseiKIN(ip,k,iq)
           if(ioptPH.ge.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >      -dek_r(iq)*dconjg(phaserKIN(ip,k,iq))
     >      +ai*dek_i(iq)*dconjg(phaseiKIN(ip,k,iq))
           elseif(ioptPH.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >      +dek_r(iq)*phaserKIN(ip,k,iq)
     >      +ai*dek_i(iq)*phaseiKIN(ip,k,iq)
           endif
          elseif(jn.lt.0) then
           jn=-jn
           matrix(ip,jn)=matrix(ip,jn)
     >     -dek_r(iq)*phaserKIN(ip,k,iq)
     >     -ai*dek_i(iq)*phaseiKIN(ip,k,iq)
           if(ioptPH.ge.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >      +dek_r(iq)*dconjg(phaserKIN(ip,k,iq))
     >      -ai*dek_i(iq)*dconjg(phaseiKIN(ip,k,iq))
           elseif(ioptPH.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >      -dek_r(iq)*phaserKIN(ip,k,iq)
     >      -ai*dek_i(iq)*phaseiKIN(ip,k,iq)
           endif
          endif
         endif
        enddo
       enddo
      enddo

      if(iesonsite.ge.1) then
       do ip=1,nh
        matrix(ip,ip+nh)=matrix(ip,ip+nh)-donsite_r(1)+ai*donsite_i(1)
        matrix(ip+nh,ip)=matrix(ip+nh,ip)-donsite_r(1)-ai*donsite_i(1)
       enddo
      endif

      do iq=1,iesBCS
       do k=1,imultis(iq)
        do ip=1,nh
         jn=ivicns(ip,k,iq)
         if(jn.gt.0) then 
          matrix(ip+nh,jn)=matrix(ip+nh,jn)
     >    +dBCS_r(iq)*phaserBCS(ip,k,iq)
     >    +ai*dBCS_i(iq)*phaseiBCS(ip,k,iq)
          matrix(ip,jn+nh)=matrix(ip,jn+nh)
     >    +dBCS_r(iq)*dconjg(phaserBCS(ip,k,iq))
     >    -ai*dBCS_i(iq)*dconjg(phaseiBCS(ip,k,iq))
         elseif(jn.lt.0) then 
          jn=-jn
          matrix(ip+nh,jn)=matrix(ip+nh,jn)
     >    -dBCS_r(iq)*phaserBCS(ip,k,iq)
     >    -ai*dBCS_i(iq)*phaseiBCS(ip,k,iq)
          matrix(ip,jn+nh)=matrix(ip,jn+nh)
     >    -dBCS_r(iq)*dconjg(phaserBCS(ip,k,iq))
     >    +ai*dBCS_i(iq)*dconjg(phaseiBCS(ip,k,iq))
         endif
        enddo
       enddo
      enddo

      if(iesaf.eq.1) then
       do ip=1,nh
        if(ioptPH.lt.0) then
         matrix(ip,ip+nh)=matrix(ip,ip+nh)+dAF(1)*alphai(ip)
         matrix(ip+nh,ip)=matrix(ip+nh,ip)+dAF(1)*dconjg(alphai(ip))
        endif
       enddo
      endif

c check 
      do ip=1,2*nh
       do iq=1,2*nh
        test1=dreal(matrix(ip,iq))-dreal(matrix(iq,ip))
        test2=dimag(matrix(ip,iq))+dimag(matrix(iq,ip))
        if(dabs(test1).gt.1.d-8.or.dabs(test2).gt.1.d-8) then
         write(6,*) 'Not Hermitian!!',test1,test2
        endif
       enddo
      enddo
      do ip=1,nh
       do iq=nh+1,2*nh
        test1=dreal(matrix(ip,iq))-dreal(matrix(iq-nh,ip+nh))
        test2=dimag(matrix(ip,iq))-dimag(matrix(iq-nh,ip+nh))
        if(dabs(test1).gt.1.d-8.or.dabs(test2).gt.1.d-8) then
         write(6,*) 'Not singlet pairing!!',test1,test2
        endif
       enddo
      enddo

c Matrix diagonalization

      ALLOCATE(psip(4*nh*nh+22*nh))
      ALLOCATE(rpsip(6*nh))

      n1=4*nh*nh+1
      call zcopy(4*nh*nh,matrix,1,psip,1)

      call zheev('V','L',2*nh,psip,2*nh,eig,psip(n1),20*nh
     >,rpsip,info)

      do i=1,2*nh
       call zcopy(2*nh,psip(2*nh*(i-1)+1),1,z(1,i),1)
      enddo

      DEALLOCATE(psip)
      DEALLOCATE(rpsip)
      DEALLOCATE(matrix)

      return
      end

