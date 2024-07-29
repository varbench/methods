      program turbohub      
      implicit none 

c *** Dimensions
      INTEGER(4) nh,nh2,nbase,ndim
      INTEGER(4) nindm,nindms,maxmulti,maxmultis
      INTEGER(4) nelup,neldo,nel,nelm2
      INTEGER(4) iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf
     >,iesslater
      INTEGER(4) np,nmat,nmats,npsov
      INTEGER(4) irange,ndimtable,ndimtot
c *** Lattice
      INTEGER(4), dimension(:), allocatable :: imulti,imultis
      INTEGER(4), dimension(:,:,:), allocatable :: ivicn,ivicns
      COMPLEX(8), dimension(:,:,:), allocatable :: phaserKIN,phaserBCS
      COMPLEX(8), dimension(:,:,:), allocatable :: phaseiKIN,phaseiBCS
c *** State
      REAL(8) muf
      REAL(8), dimension(:), allocatable :: tvar
      REAL(8), dimension(:), allocatable :: dBCS_r,dBCS_i,dek_r,dek_i
      REAL(8), dimension(:), allocatable :: dmu,donsite_r,donsite_i
      REAL(8), dimension(:), allocatable :: dAF
      REAL(8), dimension(:), allocatable :: vjspin,vjdens
      REAL(8), dimension(:), allocatable :: eig
      REAL(8), dimension(:,:), allocatable :: v,vsz
      COMPLEX(8), dimension(:), allocatable :: alphai
      COMPLEX(8), dimension(:,:), allocatable :: z
c *** Hamiltonian
      INTEGER(4) ioptPH
      REAL(8) uvv,jperp
      REAL(8), dimension(:), allocatable :: tham,jham
      REAL(8), dimension(:,:), allocatable :: thop,thopv
      REAL(8), dimension(:,:), allocatable :: jhop,jhopv
c *** Walkers
      INTEGER(2), dimension(:), allocatable :: jbraw
      INTEGER(2), dimension(:,:), allocatable :: iconf
      INTEGER(2), dimension(:,:), allocatable :: kelup,keldo
      INTEGER(4), dimension(:), allocatable :: jbra
      INTEGER(4), dimension(:), allocatable :: iconfm
      REAL(4), dimension(:), allocatable :: wconfw
      REAL(8), dimension(:), allocatable :: wconfn,wsto
      REAL(8), dimension(:), allocatable :: diag,diagfn
      REAL(8), dimension(:), allocatable :: diagpart
      REAL(8), dimension(:,:), allocatable :: tabpip,tabpipSz
      COMPLEX(8), dimension(:), allocatable :: enert,enert2
      COMPLEX(8), dimension(:,:,:), allocatable :: winv
      COMPLEX(8), dimension(:,:,:), allocatable :: table
c *** Optimization
      INTEGER(4) nwspin,nwdens,nwslater
      REAL(8) epsdgel,tau
      INTEGER(4), dimension (:), allocatable :: ipos
      REAL(8), dimension (:), allocatable :: alpha,scalpar
      REAL(8), dimension (:), allocatable :: etotw
      REAL(8), dimension (:), allocatable :: sovdt
      REAL(8), dimension (:,:,:), allocatable :: sov
      COMPLEX(8), dimension(:), allocatable :: econf
      COMPLEX(8), dimension(:,:), allocatable :: matrix
      COMPLEX(8), dimension (:,:,:), allocatable :: opmat
c *** Monte Carlo
      INTEGER(4) itestr,itest,iseed,iopt,iread,nw,ngen,nscra,nbra
     >,nweight,ibinit,nraw,nn,ng,ngs,iend,ngr,ngsr,nwr,nhr,nbrar,npr
     >,ireadr,iendr,nelupr,neldor,icdiff,nacc,ngg,ngn,inext,nleft
     >,iout,indrange,indvic,ntry,indout,irtest,kk,jn,istest,indspin
     >,jout,ioutup,ioutdo,joutup,joutdo,ioutr,jj,ii,indc,k,nint
      REAL(8) epst,etry,tbra,gamma,tacc,etryr,tbrar,wbra,sumdiff,time
     >,timescra,timep,lambda,veff,tleft,costw,pdiag,drand1,ttry
     >,ratio,upvham,ener,ener2,diagsum,wtot,rata
c *** Other
      INTEGER(4), dimension(:), allocatable :: ipsip
      REAL(8), dimension(:), allocatable :: zeta,psip
      REAL(8), dimension(:), allocatable :: rpsip
      COMPLEX(8), dimension(:,:), allocatable :: cpsip

      INTEGER(4) i,j,iscramax
      character*11 nomefile

      namelist /montecarlo/ itestr,iseed,iopt,iread,nw,ngen,nscra,epst
      namelist /fixednode/ etry,tbra,gamma 
      namelist /vmc/ nbra
      namelist /optimization/ nbra,nweight,ibinit,epsdgel,tau,nraw

      nomefile='randomseeds' !//char(0)

      itestr=0
      iseed=0
      iopt=1
      iread=0
      nw=1
      ngen=0
      nscra=100
      epst=1.d-10
      etry=0.d0
      tbra=1.d0
      gamma=0.d0
      nbra=1
      nweight=0
      ibinit=0
      epsdgel=0.d0
      tau=0.d0
      nraw=0

      read(5,montecarlo)
      if(itestr.eq.0) then
       write(6,*) 'itestr must be specified'
       stop
      endif
      if(iseed.eq.0) then
       write(6,*) 'iseed must be specified'
       stop
      endif
      if(ngen.eq.0) then
       write(6,*) 'ngen must be specified'
       stop
      endif
      write(6,*)
      write(6,*) 'itestr   =',itestr
      write(6,*) 'iopt     =',iopt
      write(6,*) 'iread    =',iread
      write(6,*) 'Nw       =',nw
      write(6,*) 'nscra    =',nscra
      if(itestr.eq.1) then
       read(5,fixednode) 
       write(6,*) '*********************************'
       write(6,*) 'Fixed-node run'
       write(6,*) '*********************************'
       if(tbra.eq.0.d0) then
        write(6,*) 'tbra must be specified'
        stop
       endif
       write(6,*) 'etry     =',etry
       write(6,*) 'tbra     =',tbra
       write(6,*) 'gamma    =',gamma
      elseif(itestr.eq.2) then
       read(5,vmc) 
       write(6,*) '*********************************'
       write(6,*) 'Variational run'
       write(6,*) '*********************************'
       write(6,*) 'nbra     =',nbra
      elseif(itestr.eq.-2) then
       read(5,optimization) 
       write(6,*) '*********************************'
       write(6,*) 'Optimization run'
       write(6,*) '*********************************'
       if(nweight.eq.0) then
        write(6,*) 'nweight must be specified'
        stop
       endif
       if(ibinit.eq.0) then
        write(6,*) 'ibinit must be specified'
        stop
       endif
       write(6,*) 'nbra     =',nbra
       write(6,*) 'nweight  =',nweight
       write(6,*) 'ibinit   =',ibinit
       write(6,*) 'epsdgel  =',epsdgel
       write(6,*) 'tau      =',tau
       write(6,*) 'nraw     =',nraw
      else
       write(6,*) 'Wrong itestr',itestr
       stop
      endif

      itest=abs(itestr)

c     itestr=2 VMC
c     itestr=1 FN
c     itestr=-2 Full Minimizator

c iopt=0   continue the run starting from a given conf 
c iopt=1   begin the calculation from scratch
c iopt=2   begin with old conf but restart averages  

      call rand_init(iseed)

      open(unit=10,file='fort.10',form='unformatted',status='unknown')
      open(unit=11,file='fort.11',form='unformatted',status='unknown')
      open(unit=12,file='fort.12',form='unformatted',status='unknown')
      if(itestr.eq.-2) open(unit=13,file='fort.13',form='unformatted'
     >,status='unknown')
      if(iread.eq.1) open(unit=15,file='fort.15',form='unformatted'
     >,status='unknown')

      rewind(10)
      rewind(12) 
      if(itestr.eq.-2) rewind(13)
      if(iread.eq.1) rewind(15)

c *** READING PART, FROM THE PREVIOUS CODE
      read(10) ndim,nbase,nh,nindm,nindms,maxmulti,maxmultis
      read(10) ioptPH,irange,nelup,neldo

      ALLOCATE(imulti(nindm))
      ALLOCATE(imultis(nindms))
      ALLOCATE(ivicn(nh,maxmulti,nindm))
      ALLOCATE(ivicns(nh,maxmultis,nindms))
      ALLOCATE(phaserKIN(nh,maxmultis,nindms))
      ALLOCATE(phaseiKIN(nh,maxmultis,nindms))
      ALLOCATE(phaserBCS(nh,maxmultis,nindms))
      ALLOCATE(phaseiBCS(nh,maxmultis,nindms))

      read(10) imulti,imultis,ivicn,ivicns
      read(10) phaserKIN,phaserBCS,phaseiKIN,phaseiBCS
      read(10) iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf,muf

      ALLOCATE(z(2*nh,2*nh))
      ALLOCATE(eig(2*nh))
      ALLOCATE(v(nh,nh))
      ALLOCATE(vsz(nh,nh))
      ALLOCATE(dBCS_r(iesBCS))
      ALLOCATE(dBCS_i(iesBCS))
      ALLOCATE(dek_r(ieskin))
      ALLOCATE(dek_i(ieskin))
      ALLOCATE(dmu(iesmu))
      ALLOCATE(donsite_r(iesonsite))
      ALLOCATE(donsite_i(iesonsite))
      ALLOCATE(dAF(iesaf))
      ALLOCATE(alphai(nh))
      ALLOCATE(vjspin(iesspin))
      ALLOCATE(vjdens(iesdens))

      read(10) z,eig,vjspin,vjdens,v,vsz,dBCS_r,dBCS_i,dek_r,dek_i,dmu
     >,donsite_r,donsite_i,dAF,alphai

      ALLOCATE(tham(irange))
      ALLOCATE(jham(irange))
      ALLOCATE(tvar(irange))

      read(10) tham,jham,tvar,jperp,uvv
c *** END READING PART, FROM THE PREVIOUS CODE

      nel=nelup+neldo
      nh2=2*nh
      nelm2=2*nh*nel

      np=iesspin+iesdens+2*iesBCS+2*ieskin+iesmu+2*iesonsite+iesaf
c the factor 2 is due to the fact that there are real and imaginary parts
      nmat=np+1
      npsov=nmat*nmat
      nmats=2*iesBCS+2*ieskin+iesmu+2*iesonsite+iesaf

      iscramax=8*nh*nh+22*nh+np*np+13*np
      ALLOCATE(ipsip(nh+2*nw+2*np))
      ALLOCATE(psip(iscramax))
      ALLOCATE(cpsip(2*nh,2*nh))
      ALLOCATE(rpsip(6*nh))
      ALLOCATE(zeta(nw+1))

      if(itestr.eq.-2) then
       ALLOCATE(alpha(nmat))
       ALLOCATE(scalpar(nmat))
       ALLOCATE(ipos(nmat))
       do i=1,nmat
        scalpar(i)=1.d0
        alpha(i)=0.d0
       enddo
       do i=1,nraw
        read(5,*) nn,tacc,(ipos(j),j=1,nn)
        do j=1,nn
         if(ipos(j).gt.np) then
          write(6,*) 'Wrong position for the parameter',ipos(j),np
          stop
         endif
         scalpar(ipos(j))=tacc
        enddo
       enddo
       if(nraw.ne.0) then
        write(6,*)
        write(6,*) 'Scalpar'
        do i=1,np
         write(6,*) i,scalpar(i)
        enddo
        write(6,*)
       endif
      endif

      nwspin=0
      nwdens=iesspin*nw
      nwslater=(iesspin+iesdens)*nw

      if(iesBCS+ieskin+iesmu+iesonsite.ne.0) then 
       iesslater=1
      else
       iesslater=0
      endif

      if(itestr.eq.-2) then
       ALLOCATE(matrix(2*nh,2*nh))
       if(iesslater.eq.1) then
        ALLOCATE(opmat(2*nh,2*nh,nmats))
        call makesovop(nh,nel,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,iesBCS,ieskin,iesmu,iesonsite,iesaf,matrix,ioptPH
     >,z,eig,cpsip,opmat)
       endif

       ALLOCATE(econf(nw*nmat))
       ALLOCATE(etotw(nmat))
       ALLOCATE(sovdt(nmat*nmat))
       ALLOCATE(sov(nmat,nmat,2))
       do j=1,nw*nmat
        econf(j)=(1.d0,0.d0)
       enddo
       do j=1,nmat
        etotw(j)=0.d0
       enddo
       do j=1,npsov
        sovdt(j)=0.d0 
       enddo
      endif

      ALLOCATE(thop(nh,nh))
      ALLOCATE(thopv(nh,nh))
      ALLOCATE(jhop(nh,nh))
      ALLOCATE(jhopv(nh,nh))

      call uphop(nh,irange,ivicn,imulti,maxmulti
     >,v,vsz,tham,jham,thop,thopv,jhop,jhopv,ioptPH)

      ALLOCATE(wconfn(nw))
      ALLOCATE(wsto(nw))
      ALLOCATE(diag(nw))
      ALLOCATE(diagfn(nw))
      ALLOCATE(enert(nw))
      ALLOCATE(enert2(nw))
      ALLOCATE(diagpart(nw))
      ALLOCATE(iconf(nh,nw))
      ALLOCATE(kelup(nh,nw))
      ALLOCATE(keldo(nh,nw))

      ndimtable=0
      do i=1,irange
       do j=1,imulti(i)
        ndimtable=ndimtable+1
       enddo
      enddo
      ndimtot=nh*ndimtable
      ALLOCATE(table(nh,ndimtable,nw))
      ALLOCATE(winv(2*nh,nel,nw))
      ALLOCATE(tabpip(nh,nw))
      ALLOCATE(tabpipSz(nh,nw))

      call dscal(2*ndimtot*nw,0.d0,table,1)
      call dscal(2*nelm2*nw,0.d0,winv,1)

      if(iopt.eq.1) then

       call initspin(nw,nh,nelup,neldo,nel,iconf,kelup,keldo
     >,wconfn,z,ipsip,cpsip,epst,ioptPH)

       ng=0
       ngs=0
       iend=0

      elseif(iopt.eq.0.or.iopt.eq.2) then

       call read_seeds(nomefile)
       rewind(11)
       read(11) ngr,ngsr,etryr,nwr,nhr,nbrar,tbrar,npr,ireadr,iendr
     >,nelupr,neldor

       if(iopt.eq.0) then 
        iread=ireadr
        etry=etryr
        tbra=tbrar
        nbra=nbrar
        iend=iendr
        ng=ngr
        ngs=ngsr
       elseif(iopt.eq.2) then
        iend=0
        ng=0
        ngs=0
       endif

       if(nhr.ne.nh) then
        write(6,*) 'The calculation must continue with L=',nhr
        stop
       endif 
       if(npr.ne.np) then
        write(6,*) 'The calculation must continue with np=',npr
        stop
       endif 
       if(nwr.ne.nw) then
        write(6,*) 'The calculation must continue with nw=',nwr
        stop
       endif
       if(nelupr.ne.nelup) then
        write(6,*) 'The calculation must continue with nelup=',nelupr
        stop
       endif
       if(neldor.ne.neldo) then
        write(6,*) 'The calculation must continue with neldo=',neldor
        stop
       endif

       read(11) ((iconf(i,j),i=1,nh),j=1,nw),(wconfn(j),j=1,nw)

       do i=1,nw
        wconfn(i)=1.d0
        call upkel(nh,nelup,kelup(1,i),keldo(1,i),iconf(1,i))
       enddo

       if(iopt.eq.0) then 
        if(itestr.eq.-2) then
         do i=1,ngs
          read(13)
         enddo
        endif
        do i=1,ng
         if(iread.eq.1) then
          read(12)
          read(15)
         else
          read(12)
         endif
        enddo
       endif

      endif

      ALLOCATE(jbraw(nw))
      ALLOCATE(jbra(nw))
      do j=1,nw
       jbraw(j)=j
       jbra(j)=j
      enddo
      if(iread.eq.1) then
       ALLOCATE(iconfm((nw*nh)/16+1))
       ALLOCATE(wconfw(nw))
      endif

      icdiff=0
      wbra=1.d0

      nacc=0
      sumdiff=0.d0

      time=mclock()/100.d0
      timescra=0.d0 

      lambda=-nh*etry

      ngg=ng   ! contatore file 12
      ngn=ngs  ! contatore file 13
      write(6,*) 'Starting ngg, ngn=',ngg,ngn

c======================================================================
c======================================================================

      i=iend
      inext=iend+nweight

      DO WHILE(i.lt.ngen+iend)

       i=i+1
       if(mod(i,nscra).eq.1.or.nscra.eq.1.or.i.eq.iend+1) then
        timep=mclock()/100.d0

        call upscratchub(nh,nw,nelup,neldo,nel
     >,tabpip,tabpipSz,diag,v,vSz,jham,jperp,uvv,irange,ivicn
     >,maxmulti,imulti,iconf,kelup,keldo,winv,ipsip,cpsip,z,ioptPH)

        if(itest.ne.2.and.i.eq.iend+1) then 
         do j=1,nw
          call uptablehubf(nh,irange,ivicn,maxmulti,imulti
     >,table(1,1,j),tabpip(1,j),tabpipSz(1,j),thop,jhop
     >,winv(1,1,j),iconf(1,j),kelup(1,j),keldo(1,j),ioptPH)

          call energy(ndimtot,table(1,1,j),diag(j),wsto(j),veff
     >,enert(j))
          diagfn(j)=diag(j)+(1.d0+gamma)*veff
         enddo
        endif 

        timescra=timescra+(mclock()/100.d0-timep)

       endif         !fine if sullo scratch

c --------------DO sui walkers--------------------------

       do j=1,nw    

        nleft=nbra
        tleft=tbra

        do while(tleft.gt.0.d0.and.nleft.gt.0) 

         if(itest.eq.1) then 
          costw=wsto(j)-lambda
          pdiag=diagfn(j)-wsto(j)
          zeta(1)=drand1()
          ttry=dlog(1.d0-zeta(1))/pdiag
          if(ttry.gt.tleft) ttry=tleft
         else
          iout=drand1()*nh+1
          indrange=irange*drand1()+1
          indvic=imulti(indrange)*drand1()+1
 
          call ratiovmc(iout,indrange,indvic,nh,ivicn,maxmulti
     >,ratio,tabpip(1,j),tabpipSz(1,j),thopv,jhopv,winv(1,1,j)
     >,iconf(1,j),kelup(1,j),keldo(1,j),ioptPH)

          zeta(1)=1.d0-drand1()
          if(ratio**2.lt.zeta(1)) then
           ntry=1
          else
           ntry=0
          endif
         endif

         if((itest.eq.1.and.ttry.ge.tleft).or.
     1      (itest.eq.2.and.ntry.eq.1)) then

          if(itest.eq.1) then
           if(ttry.ne.tleft) then 
            write(6,*) ' error impossible condition
     1                  ttry > tleft ',ttry,tleft
            stop
           endif
           if(wconfn(j).ne.0.d0) then 
            costw=dexp(costw*tleft)
            wconfn(j)=wconfn(j)*costw
           endif
           tleft=0.d0
          endif 
          if(itest.ne.1) nleft=nleft-1

         else   ! accept a new move 

          if(itest.eq.1) then
           costw=dexp(costw*ttry)
           wconfn(j)=wconfn(j)*costw
           tleft=tleft-ttry
          endif 
          if(itest.ne.1) nleft=nleft-1
          nacc=nacc+1
          if(itest.eq.1) then 
           zeta(1)=drand1()*(wsto(j)-diagfn(j))

           call random(nh,ndimtot,zeta,table(1,1,j),iout,indout,gamma)

           irtest=0
           kk=0
           do while(irtest.lt.indout)
            kk=kk+1
            irtest=irtest+imulti(kk)
           enddo

           indrange=kk
           irtest=irtest-imulti(kk)
           indvic=indout-irtest
          endif

          jn=abs(ivicn(iout,indvic,indrange))
          istest=iconf(iout,j)+iconf(jn,j)
          if(ioptPH.ge.0) then
           if(istest.eq.1) then
            indspin=1
           elseif(istest.eq.-1) then
            indspin=-1
           elseif(istest.eq.2) then
            indspin=2
           endif
          else
           if(istest.eq.1) then
            indspin=1
           elseif(istest.eq.-1) then
            indspin=-1
           elseif(istest.eq.0.and.iconf(iout,j).ne.0) then
            indspin=2
           endif
          endif

          if(iconf(iout,j).eq.0.or.iconf(jn,j).eq.0) then
           if(iconf(iout,j).eq.0) then
            jout=jn
           else
            jout=iout
            iout=jn
           endif
          elseif(iconf(iout,j).ne.0.and.iconf(jn,j).ne.0) then
           if(iconf(iout,j).lt.iconf(jn,j)) then
            jout=jn
           else
            jout=iout
            iout=jn
           endif
          endif

          if(indspin.eq.2) then
           if(ioptPH.ge.0) then
            ioutup=iout
            ioutdo=iout
            joutup=jout
            joutdo=jout
           else
            joutup=jout
            ioutup=iout
            joutdo=iout
            ioutdo=jout
           endif
          endif

          IF(indspin.eq.1) then
       call upinvhop(nh,nel,iout,jout,kelup(jout,j),winv(1,1,j)
     >,cpsip,epst)
       call uptabtothubf(nh,nel,tabpip(1,j),tabpipSz(1,j)
     >,iout,jout,iconf(1,j),v(1,iout),v(1,jout),vSz(1,iout),vSz(1,jout)
     >,kelup(1,j),1,ioptPH)
          ELSEIF(indspin.eq.-1) then
       call upinvhop(nh,nel,iout+nh,jout+nh,keldo(jout,j),winv(1,1,j)
     >,cpsip,epst)
       call uptabtothubf(nh,nel,tabpip(1,j),tabpipSz(1,j)
     >,iout,jout,iconf(1,j),v(1,iout),v(1,jout),vSz(1,iout),vSz(1,jout)
     >,keldo(1,j),-1,ioptPH)
          ELSEIF(indspin.eq.2) then
       ioutr=ioutdo+nh
       call upinvhop2(nh2,nel,ioutup,kelup(joutup,j)
     >,ioutr,keldo(joutdo,j),winv(1,1,j),cpsip,epst)
       call uptabtothubf(nh,nel,tabpip(1,j),tabpipSz(1,j)
     >,ioutup,joutup,iconf(1,j),v(1,ioutup),v(1,joutup),vSz(1,ioutup)
     >,vSz(1,joutup),kelup(1,j),1,ioptPH)
       call uptabtothubf(nh,nel,tabpip(1,j),tabpipSz(1,j)
     >,ioutdo,joutdo,iconf(1,j),v(1,ioutdo),v(1,joutdo),vSz(1,ioutdo)
     >,vSz(1,joutdo),keldo(1,j),-1,ioptPH)
          ENDIF

          diag(j)=0.d0
          if(ioptPH.ge.0) then
           do jj=1,irange
            do kk=1,imulti(jj)
             do ii=1,nh
              jn=abs(ivicn(ii,kk,jj))
              if(jn.ne.0) then
               if(iconf(ii,j)+iconf(jn,j).eq.2)
     >         diag(j)=diag(j)+0.5d0*jham(jj)*jperp
              endif
             enddo
            enddo
           enddo
          else
           do jj=1,irange
            do kk=1,imulti(jj)
             do ii=1,nh
              jn=abs(ivicn(ii,kk,jj))
              if(jn.ne.0) then
               if(iconf(ii,j)*iconf(jn,j).eq.-1)
     >         diag(j)=diag(j)+0.5d0*jham(jj)*jperp
              endif
             enddo
            enddo
           enddo
          endif
          if(uvv.ne.0) diag(j)=diag(j)-
     >    upvham(uvv,jperp,jham,nh,irange,ivicn,maxmulti,imulti
     >,iconf(1,j),ioptPH)

          diag(j)=0.5d0*diag(j)

          if(itest.eq.1) then 
           call uptablehubf(nh,irange,ivicn,maxmulti,imulti
     >,table(1,1,j),tabpip(1,j),tabpipSz(1,j),thop,jhop
     >,winv(1,1,j),iconf(1,j),kelup(1,j),keldo(1,j),ioptPH)

           call energy(ndimtot,table(1,1,j),diag(j),wsto(j),veff
     >,enert(j))
           diagfn(j)=diag(j)+(1.d0+gamma)*veff
          endif

         endif    ! end if accepted meas

        enddo ! end do for the tleft or nleft

       enddo ! end do for the walkers

c Calculation of observables
       do j=1,nw
        call uptablehubf(nh,irange,ivicn,maxmulti,imulti
     >,table(1,1,j),tabpip(1,j),tabpipSz(1,j),thop,jhop
     >,winv(1,1,j),iconf(1,j),kelup(1,j),keldo(1,j),ioptPH)

        IF(itestr.eq.-2) THEN
         if(iesspin.ne.0) then
          call upvpotz(iesspin,nw,nh,ivicn,maxmulti,imulti,iconf(1,j)
     >,econf(nwspin+j),ioptPH)
         endif
         if(iesdens.ne.0) then 
          call upvpot(iesdens,nw,nh,ivicn,maxmulti,imulti,iconf(1,j)
     >,econf(nwdens+j),ioptPH)
         endif

         if(iesslater.eq.1) then
          call upallcorr(nw,nh2,nmats,opmat,winv(1,1,j)
     >,kelup(1,j),keldo(1,j),econf(nwslater+j))
         endif

         econf(np*nw+j)=(1.d0,0.d0)
        ENDIF

        call energy(ndimtot,table(1,1,j),diag(j),wsto(j),veff
     >,enert(j))
        enert(j)=enert(j)/dble(nh)

        enert2(j)=enert(j)*dconjg(enert(j))
        diagpart(j)=-diag(j)/dble(nh)
       enddo

       ener=enert(1)*wconfn(1)
       ener2=enert2(1)*wconfn(1)
       diagsum=diagpart(1)*wconfn(1)
       wtot=wconfn(1)
       do j=2,nw
        wtot=wtot+wconfn(j)
        ener=ener+enert(j)*wconfn(j)
        ener2=ener2+enert2(j)*wconfn(j)
        diagsum=diagsum+diagpart(j)*wconfn(j)
       enddo
       ener=ener/wtot
       ener2=ener2/wtot
       diagsum=diagsum/wtot

c---------------------------------------------------------------------

       if(iread.eq.1) then
        do kk=1,nw
         wconfw(kk)=wconfn(kk)
        enddo
       endif

       if(itest.eq.1) then 

        zeta(nw+1)=drand1()
        call branchingo(nw,wconfn,wbra,zeta,icdiff,ipsip,jbra)
        sumdiff=sumdiff+icdiff

        do j=1,nw
         jbraw(j)=jbra(j)
        enddo

       elseif(itestr.eq.-2) then 

        IF(i.eq.inext) THEN 

         call reweight0(nw,np,nmat,etotw,ipsip,psip,alpha,sov
     >,sovdt,econf,epsdgel,1,nh,enert,scalpar)

         indc=0
         if(iesspin.ne.0) then 
          do k=indc+1,indc+iesspin
           vjspin(k-indc)=vjspin(k-indc)+alpha(k)*tau
           alpha(k)=vjspin(k-indc)
          enddo
         endif
         indc=indc+iesspin
         if(iesdens.ne.0) then 
          do k=indc+1,indc+iesdens
           vjdens(k-indc)=vjdens(k-indc)+alpha(k)*tau
           alpha(k)=vjdens(k-indc)
          enddo
         endif
         indc=indc+iesdens
         if(iesBCS.ne.0) then 
          do k=indc+1,indc+iesBCS
           dBCS_r(k-indc)=dBCS_r(k-indc)+alpha(k)*tau
           alpha(k)=dBCS_r(k-indc)
          enddo
          indc=indc+iesBCS
          do k=indc+1,indc+iesBCS
           dBCS_i(k-indc)=dBCS_i(k-indc)+alpha(k)*tau
           alpha(k)=dBCS_i(k-indc)
          enddo
         endif
         indc=indc+iesBCS
         if(ieskin.ne.0) then 
          do k=indc+1,indc+ieskin
           dek_r(k-indc)=dek_r(k-indc)+alpha(k)*tau
           alpha(k)=dek_r(k-indc)
          enddo
          indc=indc+ieskin
          do k=indc+1,indc+ieskin
           dek_i(k-indc)=dek_i(k-indc)+alpha(k)*tau
           alpha(k)=dek_i(k-indc)
          enddo
         endif
         indc=indc+ieskin
         if(iesmu.ne.0) then 
          do k=indc+1,indc+iesmu
           dmu(k-indc)=dmu(k-indc)+alpha(k)*tau
           alpha(k)=dmu(k-indc)
          enddo
         endif
         indc=indc+iesmu
         if(iesonsite.ne.0) then 
          do k=indc+1,indc+iesonsite
           donsite_r(k-indc)=donsite_r(k-indc)+alpha(k)*tau
           alpha(k)=donsite_r(k-indc)
          enddo
         endif
         indc=indc+iesonsite
         if(iesonsite.ne.0) then 
          do k=indc+1,indc+iesonsite
           donsite_i(k-indc)=donsite_i(k-indc)+alpha(k)*tau
           alpha(k)=donsite_i(k-indc)
          enddo
         endif
         indc=indc+iesonsite
         if(iesaf.ne.0) then 
          do k=indc+1,indc+iesaf
           dAF(k-indc)=dAF(k-indc)+alpha(k)*tau
           alpha(k)=dAF(k-indc)
          enddo
         endif

         call ekdk(nh,irange,ivicn,ivicns,maxmulti,maxmultis
     >,imulti,imultis,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS
     >,iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf,iesslater
     >,vjspin,vjdens,dBCS_r,dBCS_i,dek_r,dek_i,dmu,donsite_r,donsite_i
     >,dAF,alphai,cpsip,rpsip,matrix,vsz,v,z,muf,tvar,eig,ioptPH)

         if(iesslater.eq.1)
     > call makesovop(nh,nel,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,iesBCS,ieskin,iesmu,iesonsite,iesaf,matrix,ioptPH
     >,z,eig,cpsip,opmat)

         call uphop(nh,irange,ivicn,imulti,maxmulti
     >,v,vsz,tham,jham,thop,thopv,jhop,jhopv,ioptPH)

         call upscratchub(nh,nw,nelup,neldo,nel
     >,tabpip,tabpipSz,diag,v,vSz,jham,jperp,uvv,irange,ivicn
     >,maxmulti,imulti,iconf,kelup,keldo,winv,ipsip,cpsip,z,ioptPH)

         write(6,*) ' par   =',i,(alpha(k),k=1,nmat)

         ngn=ngn+1
         write(13) (alpha(k),k=1,nmat)

         do k=1,nmat
          etotw(k)=0.d0
         enddo
         do k=1,npsov
          sovdt(k)=0.d0
         enddo

         inext=inext+nweight

        ELSE  ! i.eq.inext

         call reweight0(nw,np,nmat,etotw,ipsip,psip,alpha,sov
     >,sovdt,econf,epsdgel,0,nh,enert,scalpar)

         if(i.eq.inext+ibinit-nweight) then 
          do k=1,nmat
           etotw(k)=0.d0
          enddo
          do k=1,npsov
           sovdt(k)=0.d0
          enddo
         endif

        ENDIF

       endif

       if(itest.eq.1) then
        if(iread.eq.1) then 
         call convert(nh,nw,nint,iconf,iconfm,ioptPH)
        endif
        call reshuffhub(nh,ndimtot,nw,jbra,iconf
     >,table,tabpip,tabpipSz,winv,nelm2,wsto
     >,diagfn,diag,kelup,keldo,enert)
       endif
       if(itest.eq.2) then
        wbra=wtot
        if(iread.eq.1) then 
         call convert(nh,nw,nint,iconf,iconfm,ioptPH)
        endif
       endif

       if(iread.eq.1) then 
        ngg=ngg+1
        write(12) wbra,wtot,ener,ener2,diagsum
        write(15) wbra,(iconfm(k),k=1,nint)
     >,(jbraw(k),k=1,nw),(wconfw(k),k=1,nw)
       else
        ngg=ngg+1
        write(12) wbra,wtot,ener,ener2,diagsum
       endif
       
      enddo ! enddo for the generations

      iend=i

      if(itestr.eq.-2) then 

       write(6,*) ' Write the parameters of the final wave function '

       rewind(10)
       write(10) ndim,nbase,nh,nindm,nindms,maxmulti,maxmultis
       write(10) ioptPH,irange,nelup,neldo
       write(10) imulti,imultis,ivicn,ivicns
       write(10) phaserKIN,phaserBCS,phaseiKIN,phaseiBCS
       write(10) iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf,muf
       write(10) z,eig,vjspin,vjdens,v,vsz,dBCS_r,dBCS_i,dek_r,dek_i,dmu
     >,donsite_r,donsite_i,dAF,alphai
       write(10) tham,jham,tvar,jperp,uvv

      endif

      call write_seeds(nomefile)

      rewind(11)
      write(11) ngg,ngn,etry,nw,nh,nbra,tbra,np,iread,iend
     >,nelup,neldo

      write(11) ((iconf(k,j),k=1,nh),j=1,nw),(wconfn(j),j=1,nw)

      close(10)
      close(11) 
      close(12)
      if(itestr.eq.-2) close(13)
      if(iread.eq.1) close(15)

      time=mclock()/100.d0-time
      sumdiff=sumdiff/dble(nw*ngen)
      rata=dble(nacc)/(ngen*nbra*tbra*nw)

      write(6,*) ' Total time =',time
      time=time-timescra
      timescra=timescra/time
      write(6,*) ' accept. rate  off diagonal moves =',rata
      write(6,*) '# of ind. walkers/ # walkers ',sumdiff
      write(6,*) ' Ratio time for scratch update =',timescra

      DEALLOCATE(imulti)
      DEALLOCATE(imultis)
      DEALLOCATE(ivicn)
      DEALLOCATE(ivicns)
      DEALLOCATE(phaserKIN)
      DEALLOCATE(phaseiKIN)
      DEALLOCATE(phaserBCS)
      DEALLOCATE(phaseiBCS)
      DEALLOCATE(z)
      DEALLOCATE(eig)
      DEALLOCATE(v)
      DEALLOCATE(vsz)
      DEALLOCATE(dBCS_r)
      DEALLOCATE(dBCS_i)
      DEALLOCATE(dek_r)
      DEALLOCATE(dek_i)
      DEALLOCATE(dmu)
      DEALLOCATE(donsite_r)
      DEALLOCATE(donsite_i)
      DEALLOCATE(dAF)
      DEALLOCATE(alphai)
      DEALLOCATE(vjspin)
      DEALLOCATE(vjdens)
      DEALLOCATE(tham)
      DEALLOCATE(jham)
      DEALLOCATE(tvar)
      DEALLOCATE(ipsip)
      DEALLOCATE(psip)
      DEALLOCATE(cpsip)
      DEALLOCATE(rpsip)
      DEALLOCATE(zeta)
      DEALLOCATE(thop)
      DEALLOCATE(thopv)
      DEALLOCATE(jhop)
      DEALLOCATE(jhopv)
      DEALLOCATE(wconfn)
      DEALLOCATE(wsto)
      DEALLOCATE(diag)
      DEALLOCATE(diagfn)
      DEALLOCATE(enert)
      DEALLOCATE(enert2)
      DEALLOCATE(diagpart)
      DEALLOCATE(iconf)
      DEALLOCATE(kelup)
      DEALLOCATE(keldo)
      DEALLOCATE(table)
      DEALLOCATE(winv)
      DEALLOCATE(tabpip)
      DEALLOCATE(tabpipSz)
      DEALLOCATE(jbraw)
      DEALLOCATE(jbra)

      if(iread.eq.1) then
       DEALLOCATE(iconfm)
       DEALLOCATE(wconfw)
      endif

      if(itestr.eq.-2) then
       DEALLOCATE(matrix)
       if(iesslater.eq.1) then
        DEALLOCATE(opmat)
       endif
       DEALLOCATE(alpha)
       DEALLOCATE(scalpar)
       DEALLOCATE(ipos)
       DEALLOCATE(econf)
       DEALLOCATE(etotw)
       DEALLOCATE(sovdt)
       DEALLOCATE(sov)
      endif

      stop 
      end

c List of subroutines
c
c convert
c random
c energy
c initspin
c upkel
c upscratchub
c uptabpip
c uptabtothubf
c upinvhop
c ratiovmc
c uptablehubf
c
c branchingo
c upjbra
c reshuffhub
c
c reweight0
c dgelssn
c mscale
c dscalzero
c makesovop
c upmatrix
c upallcorr
c upvpot
c upvpotz
c uphop
c ekdk
c upvham
c upinvhop2
c     ------------------------------------------------------
c     new convert
c     ------------------------------------------------------

      subroutine convert(L,nw,nwords,iconf,iaux,iopt)
      implicit none
      integer*4 L,nw,iopt
      integer*2 iconf(*)
      integer*4 iaux(*)
      integer*4 ndim,nwords,nrest,i,j,imax,ind

      ndim=L*nw
      nwords=ndim/16
      nrest=ndim-16*nwords
      
      if(nrest.ne.0) nwords=nwords+1
      
      do i=1,nwords
       iaux(i)=0
      enddo

      do j=1,nwords
       if(j.ne.nwords.or.nrest.eq.0) then
        imax=16
       else
        imax=nrest
       endif
       if(iopt.ge.0) then 
        do i=0,imax-1
         ind=16*(j-1)+i+1
         if(iconf(ind).eq.2) iaux(j)=ibset(iaux(j),2*i)
         if(iconf(ind).eq.0) iaux(j)=ibset(iaux(j),2*i+1)
         if(iconf(ind).eq.1) then
          iaux(j)=ibset(iaux(j),2*i)
          iaux(j)=ibset(iaux(j),2*i+1)
         endif
        enddo
       else
        do i=0,imax-1
         ind=16*(j-1)+i+1
         if(iconf(ind).eq.1) iaux(j)=ibset(iaux(j),2*i)
         if(iconf(ind).eq.-1) iaux(j)=ibset(iaux(j),2*i+1)
         if(iconf(ind).eq.2) then
          iaux(j)=ibset(iaux(j),2*i)
          iaux(j)=ibset(iaux(j),2*i+1)
         endif
        enddo
       endif
      enddo

      return
      end

      subroutine random(L,Lz,ztry,table,iout,indout,gamma)
      implicit none
      integer*4 i,L,Lz,iout,indout,indz
      real*8 try,ztry,gamma,cost
      complex(8) table(*)

      try=0.d0
      i=0
      dowhile(ztry.ge.try.and.i.lt.Lz)
       i=i+1
       if(dreal(table(i)).gt.0.d0) then
        try=try+table(i)
        indz=i
       else
        cost=-gamma*table(i)
        try=try+cost
        if(cost.ne.0.d0) indz=i
       endif
      enddo

      if(i.eq.Lz) i=indz 
      indout=(i-1)/L+1
      iout=i-(indout-1)*L

      return
      end

      subroutine energy(Lz,table,diag,wtot,veff,ener)
      implicit none
      integer*4 Lz,i
      real*8 diag,wtot,veff
      complex(8) table(*),ener

      wtot=diag
      ener=-dcmplx(diag,0.d0)
      veff=0.d0
      do i=1,Lz
       wtot=wtot+table(i)
       ener=ener-table(i)
       if(dreal(table(i)).lt.0.d0) veff=veff+table(i)
      enddo

      return
      end

      subroutine initspin(nw,L,nelup,neldo,nel,iconf,kelup,keldo
     >,wconf,wpot,ipsip,psip,epst,ioptPH)
      implicit none
      integer*4 nw,L,nelup,neldo,nel,i,j,kk,ih,pos,ioptPH
      integer*2 iconf(L,*),kelup(L,*),keldo(L,*),keli
      integer*4 ipsip(L),info,igen,iref
      real*8 wconf(*),drand1
      real*8 det(2),cdet,cost,epst,dabs
      complex(8) wpot(2*L,*),psip(L,L)

      igen=0
      iref=0
      do kk=1,nw
       cdet=0.d0
11     continue
       do i=1,L
        iconf(i,kk)=0
       enddo
       ih=0
       wconf(kk)=1.d0

       if(ioptPH.ge.0) then
        do while(ih.lt.neldo)
         pos=drand1()*L+1
         if(iconf(pos,kk).eq.0) then
          iconf(pos,kk)=-1
          ih=ih+1
         endif
        enddo
        ih=0
        do while(ih.lt.nelup)
         pos=drand1()*L+1
         if(iconf(pos,kk).eq.-1) then
          iconf(pos,kk)=2
          ih=ih+1
         endif
        enddo
       else
        do while(ih.lt.neldo)
         pos=drand1()*L+1
         if(iconf(pos,kk).eq.0) then
          iconf(pos,kk)=-1
          ih=ih+1
         endif
        enddo
        ih=0
        do while(ih.lt.nelup)
         pos=drand1()*L+1
         if(iconf(pos,kk).eq.0) then
          iconf(pos,kk)=1
          ih=ih+1
         endif
        enddo
       endif

       call upkel(L,nelup,kelup(1,kk),keldo(1,kk),iconf(1,kk))

       do i=1,L
        if(iconf(i,kk).eq.1.or.iconf(i,kk).eq.2) then
         keli=kelup(i,kk)
         do j=1,nel
          psip(keli,j)=wpot(i,j)
         enddo
        endif
        if(iconf(i,kk).eq.-1.or.iconf(i,kk).eq.2) then
         keli=keldo(i,kk)
         do j=1,nel
          psip(keli,j)=wpot(i+L,j)
         enddo
        endif
       enddo

       igen=igen+1
       call zgetrf(nel,nel,psip,L,ipsip,info)
       if(info.ne.0) then
        iref=iref+1
        goto 11
       endif
       cdet=1.d0
       do i=1,nel
        cost=abs(psip(i,i))
        if(cost.lt.cdet) cdet=cost
       enddo
       if(cdet.lt.epst) then
c       write(6,*) 'refused singular =',cdet
        iref=iref+1
        goto 11
       endif
      enddo

      if(iref.eq.0) then
       write(6,*) 'all the configurations accepted'
      else
       write(6,*) '% refused configurations = ',dble(iref)/dble(igen)
      endif

      return
      end

      subroutine upkel(L,nelup,kelup,keldo,iconf)
      implicit none
      integer*4 L,docnt,upcnt,i,nelup
      integer*2 kelup(L),keldo(L),iconf(L)

      upcnt=0
      docnt=0
      do i=1,L
       if(iconf(i).eq.1) then
        upcnt=upcnt+1
        kelup(i)=upcnt
        keldo(i)=0
       elseif(iconf(i).eq.-1) then
        docnt=docnt+1
        keldo(i)=docnt
        kelup(i)=0
       elseif(iconf(i).eq.2) then
        upcnt=upcnt+1
        docnt=docnt+1
        kelup(i)=upcnt
        keldo(i)=docnt
       else
        kelup(i)=0
        keldo(i)=0
       endif
       if(keldo(i).ne.0) keldo(i)=keldo(i)+nelup
      enddo

      return
      end

      subroutine upscratchub(L,nw,nelup,neldo,nel
     >,tabpip,tabpipSz,diag,v,vSz,jham,jperp,uvv,irange,ivicn
     >,maxmulti,imulti,iconf,kelup,keldo,winv,ipsip,psip,z,ioptPH)
      implicit none
      integer*4 nw,L,nelup,neldo,nel,i,j,k,jn,ii,info,ioptPH
      integer*4 maxmulti,irange
      integer*4 ipsip(L)
      integer*4 ivicn(L,maxmulti,*),imulti(*)
      integer*2 iconf(L,*),kelup(L,*),keldo(L,*),keli
      real*8 tabpip(L,*),tabpipSz(L,*)
     >,v(L,L),vSz(L,L),diag(*),jham(*)
     >,jperp,uvv,upvham
      complex(8) winv(2*L,nel,*),z(2*L,nel),psip(L,*)

      do k=1,nw
       if(iconf(1,k).eq.1) then
        do i=1,L
         if(ioptPH.ge.0) then 
          tabpip(i,k)=v(i,1)
          tabpipSz(i,k)=vSz(i,1)
         else
          tabpip(i,k)=v(i,1)
          tabpipSz(i,k)=vSz(i,1)
         endif
        enddo
       elseif(iconf(1,k).eq.-1) then
        do i=1,L
         if(ioptPH.ge.0) then 
          tabpip(i,k)=1.d0/v(i,1)
          tabpipSz(i,k)=vSz(i,1)
         else
          tabpip(i,k)=v(i,1)
          tabpipSz(i,k)=1.d0/vSz(i,1)
         endif
        enddo
       elseif(iconf(1,k).eq.2) then
        do i=1,L
         if(ioptPH.ge.0) then 
          tabpip(i,k)=1.d0
          tabpipSz(i,k)=vSz(i,1)**2
         else
          tabpip(i,k)=v(i,1)**2
          tabpipSz(i,k)=1.d0
         endif
        enddo
       elseif(iconf(1,k).eq.0) then
        do i=1,L
         if(ioptPH.ge.0) then 
          tabpip(i,k)=1.d0
          tabpipSz(i,k)=1.d0
         else
          tabpip(i,k)=1.d0
          tabpipSz(i,k)=1.d0
         endif
        enddo
       endif
      enddo

      do k=1,nw
       do j=2,L
        if(iconf(j,k).eq.1) then
         do i=1,L
          if(ioptPH.ge.0) then 
           tabpip(i,k)=tabpip(i,k)*v(i,j)
           tabpipSz(i,k)=tabpipSz(i,k)*vSz(i,j)
          else
           tabpip(i,k)=tabpip(i,k)*v(i,j)
           tabpipSz(i,k)=tabpipSz(i,k)*vSz(i,j)
          endif
         enddo
        elseif(iconf(j,k).eq.-1) then
         do i=1,L
          if(ioptPH.ge.0) then 
           tabpip(i,k)=tabpip(i,k)/v(i,j)
           tabpipSz(i,k)=tabpipSz(i,k)*vSz(i,j)
          else
           tabpip(i,k)=tabpip(i,k)*v(i,j)
           tabpipSz(i,k)=tabpipSz(i,k)/vSz(i,j)
          endif
         enddo
        elseif(iconf(j,k).eq.2) then 
         do i=1,L
          if(ioptPH.ge.0) then 
           tabpipSz(i,k)=tabpipSz(i,k)*vSz(i,j)**2
          else
           tabpip(i,k)=tabpip(i,k)*v(i,j)**2
          endif
         enddo
        elseif(iconf(j,k).eq.0) then 
         do i=1,L
          if(ioptPH.lt.0) then 
          endif
         enddo
        endif
       enddo
      enddo

      do k=1,nw
       do i=1,nel
        do j=1,nel
         psip(i,j)=dcmplx(0.d0,0.d0)
        enddo
       enddo

       do i=1,L
        if(iconf(i,k).eq.1.or.iconf(i,k).eq.2) then
         keli=kelup(i,k)
         do j=1,nel
          psip(keli,j)=z(i,j)
         enddo
        endif
        if(iconf(i,k).eq.-1.or.iconf(i,k).eq.2) then
         keli=keldo(i,k)
         do j=1,nel
          psip(keli,j)=z(i+L,j)
         enddo
        endif
       enddo

       call zgetrf(nel,nel,psip,L,ipsip,info)
       if(info.ne.0) then
        write(6,*) ' kelup =',(kelup(i,k),i=1,L)
        write(6,*) ' keldo =',(keldo(i,k),i=1,L)
        write(6,*) 'ERROR IN DGETRF info =',info,'nw=',k
       endif

       if(info.eq.0) then 
        call zgetri(nel,psip,L,ipsip,psip(1,nel+1),nel,info)
        if(info.ne.0) write(6,*) 'ERROR IN DGETRI info =',info
        call zgemm('N','N',2*L,nel,nel,(1.d0,0.d0),z,2*L,psip
     >,L,(0.d0,0.d0),winv(1,1,k),2*L)
       endif
      enddo

      do k=1,nw
       diag(k)=0.d0
      enddo

      do j=1,irange
       do ii=1,imulti(j)
        do i=1,L
         jn=abs(ivicn(i,ii,j))
         if(jn.ne.0) then
          if(ioptPH.ge.0) then 
           do k=1,nw
            if(iconf(i,k)+iconf(jn,k).eq.2) 
     >      diag(k)=diag(k)+0.5d0*jham(j)*jperp
           enddo
          else
           do k=1,nw
            if(iconf(i,k)*iconf(jn,k).eq.-1) 
     >      diag(k)=diag(k)+0.5d0*jham(j)*jperp
           enddo
          endif
         endif
        enddo
       enddo
      enddo
 
      if(uvv.ne.0.d0) then
       do k=1,nw
        diag(k)=diag(k)-upvham(uvv,jperp,jham,L,irange,ivicn,maxmulti
     >,imulti,iconf(1,k),ioptPH)
       enddo
      endif

      do k=1,nw
       diag(k)=0.5d0*diag(k)
      enddo

      return
      end

      subroutine uptabpip(L,dsz,tabpip,v)
      implicit none
      integer*4 L,dsz,i
      real*8 v(*),tabpip(*)

      if(dsz.eq.1) then
       do i=1,L
        tabpip(i)=tabpip(i)*v(i)
       enddo
      elseif(dsz.eq.-1) then 
       do i=1,L
        tabpip(i)=tabpip(i)/v(i)
       enddo
      endif

      return
      end

      subroutine uptabtothubf(L,nel,tabpip,tabpipSz
     >,i,j,iconf,vi,vj,vSzi,vSzj,kel,indspin,ioptPH)
      implicit none
      integer*4 L,i,j,k,szi,szj,nel,nelup,indspin,kels,ioptPH
      integer*2 iconf(*),kel(L)
      real*8 tabpip(*),tabpipSz(*)
     >,vi(*),vj(*),vSzi(*),vSzj(*)

      szi=iconf(i)
      szj=iconf(j)
      if(abs(iconf(i)).eq.1) szi=iconf(i)*indspin
      if(abs(iconf(j)).eq.1) szj=iconf(j)*indspin

      if(szi.eq.0.and.szj.eq.1) then
       if(ioptPH.ge.0) then 
        call uptabpip(L,indspin,tabpip,vi)
        call uptabpip(L,-indspin,tabpip,vj)
        call uptabpip(L,1,tabpipSz,vSzi)
        call uptabpip(L,-1,tabpipSz,vSzj)
       else
        call uptabpip(L,1,tabpip,vi)
        call uptabpip(L,-1,tabpip,vj)
        call uptabpip(L,indspin,tabpipSz,vSzi)
        call uptabpip(L,-indspin,tabpipSz,vSzj)
       endif
       kels=kel(i)
       kel(i)=kel(j)
       kel(j)=kels
       iconf(i)=iconf(j)
       iconf(j)=0
      elseif(szi.eq.0.and.szj.eq.2) then
       if(ioptPH.ge.0) then 
        call uptabpip(L,indspin,tabpip,vi)
        call uptabpip(L,-indspin,tabpip,vj)
        call uptabpip(L,1,tabpipSz,vSzi)
        call uptabpip(L,-1,tabpipSz,vSzj)
       else
        call uptabpip(L,1,tabpip,vi)
        call uptabpip(L,-1,tabpip,vj)
        call uptabpip(L,indspin,tabpipSz,vSzi)
        call uptabpip(L,-indspin,tabpipSz,vSzj)
       endif
       kels=kel(i)
       kel(i)=kel(j)
       kel(j)=kels
       iconf(i)=indspin
       iconf(j)=-iconf(i)
      elseif(szi.eq.-1.and.szj.eq.1) then
       if(ioptPH.ge.0) then
        call uptabpip(L,indspin,tabpip,vi)
        call uptabpip(L,-indspin,tabpip,vj)
        call uptabpip(L,1,tabpipSz,vSzi)
        call uptabpip(L,-1,tabpipSz,vSzj)
       else
        call uptabpip(L,1,tabpip,vi)
        call uptabpip(L,-1,tabpip,vj)
        call uptabpip(L,indspin,tabpipSz,vSzi)
        call uptabpip(L,-indspin,tabpipSz,vSzj)
       endif
       kels=kel(i)
       kel(i)=kel(j)
       kel(j)=kels
       iconf(i)=2
       iconf(j)=0
      elseif(szi.eq.-1.and.szj.eq.2) then
       if(ioptPH.ge.0) then 
        call uptabpip(L,indspin,tabpip,vi)
        call uptabpip(L,-indspin,tabpip,vj)
        call uptabpip(L,1,tabpipSz,vSzi)
        call uptabpip(L,-1,tabpipSz,vSzj)
       else
        call uptabpip(L,1,tabpip,vi)
        call uptabpip(L,-1,tabpip,vj)
        call uptabpip(L,indspin,tabpipSz,vSzi)
        call uptabpip(L,-indspin,tabpipSz,vSzj)
       endif
       kels=kel(i)
       kel(i)=kel(j)
       kel(j)=kels
       iconf(j)=iconf(i)
       iconf(i)=2
      endif

      return
      end

      subroutine upinvhop(L,nel,isite,jsite,kelj,winv,psip,epst)
      implicit none
      integer*2 kelj
      integer*4 nel,nelm,i,j,isite,jsite,L,L2,keli
      real*8 epst
      complex(8) winv(2*L,*),g,psip(2*L,*)

      L2=2*L

      g=winv(isite,kelj)

      if(abs(g).lt.epst) then 
       return
      endif
      do j=1,L2
       psip(j,1)=winv(j,kelj)/g
      enddo
      do j=1,nel
       psip(j,2)=winv(isite,j)-winv(jsite,j)
      enddo
      call zgeru(L2,nel,(-1.d0,0.d0),psip,1,psip(1,2),1,winv,L2)

      return
      end

      subroutine ratiovmc(i,indrange,indvic,L,ivicn,maxmulti
     >,table,tabpip,tabpipSz,thop,jhop,winv
     >,iconf,kelup,keldo,ioptPH)
      implicit none
      integer*4 L,i,indrange,indvic,maxmulti,jn,ioptPH,kj,ki,il,jnl
      integer*4 ivicn(L,maxmulti,*)
      integer*2 iconf(L),kelup(L),keldo(L)
      real*8 table,tabpip(*),tabpipSz(*),thop(L,L),jhop(L,L)
      complex(8) winv(2*L,*),cdet,tablec

      jn=abs(ivicn(i,indvic,indrange))

      IF(jn.eq.0) THEN
       tablec=dcmplx(0.d0,0.d0)
      ELSE
       if(ioptPH.ge.0) then 

        if(iconf(i).eq.iconf(jn)) then
         tablec=dcmplx(0.d0,0.d0)
        elseif(iconf(i).eq.0.and.iconf(jn).eq.-1) then
         kj=keldo(jn)
         il=i+L
         cdet=winv(il,kj)
         tablec=cdet*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*thop(i,jn)
        elseif(iconf(i).eq.-1.and.iconf(jn).eq.0) then
         ki=keldo(i)
         jnl=jn+L
         cdet=winv(jnl,ki)
         tablec=cdet*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*thop(i,jn)
        elseif(iconf(i).eq.-1.and.iconf(jn).eq.2) then
         kj=kelup(jn)
         cdet=winv(i,kj)
         tablec=cdet*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*thop(i,jn)
        elseif(iconf(i).eq.2.and.iconf(jn).eq.-1) then
         ki=kelup(i)
         cdet=winv(jn,ki)
         tablec=cdet*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*thop(i,jn)
        elseif(iconf(i).eq.0.and.iconf(jn).eq.2) then
         kj=kelup(jn)
         ki=keldo(jn)
         il=i+L
         cdet=winv(i,kj)*winv(il,ki)-winv(il,kj)*winv(i,ki)
         tablec=cdet*(tabpipSz(i)/tabpipSz(jn))**2*jhop(i,jn)
        elseif(iconf(i).eq.2.and.iconf(jn).eq.0) then
         ki=kelup(i)
         kj=keldo(i)
         jnl=jn+L
         cdet=winv(jn,ki)*winv(jnl,kj)-winv(jnl,ki)*winv(jn,kj)
         tablec=cdet*(tabpipSz(jn)/tabpipSz(i))**2*jhop(i,jn)
        endif

       else

        if(iconf(i).eq.iconf(jn)) then
         tablec=dcmplx(0.d0,0.d0)
        elseif(iconf(i).eq.0.and.iconf(jn).eq.1) then
         kj=kelup(jn)
         cdet=winv(i,kj)
         tablec=cdet*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*thop(i,jn)
        elseif(iconf(i).eq.1.and.iconf(jn).eq.0) then
         ki=kelup(i)
         cdet=winv(jn,ki)
         tablec=cdet*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*thop(i,jn)
        elseif(iconf(i).eq.0.and.iconf(jn).eq.-1) then
         kj=keldo(jn)
         il=i+L
         cdet=winv(il,kj)
         tablec=cdet*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*thop(i,jn)
        elseif(iconf(i).eq.-1.and.iconf(jn).eq.0) then
         ki=keldo(i)
         jnl=jn+L
         cdet=winv(jnl,ki)
         tablec=cdet*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*thop(i,jn)
        elseif(iconf(i).eq.-1.and.iconf(jn).eq.1) then
         kj=kelup(jn)
         ki=keldo(i)
         jnl=jn+L
         cdet=winv(i,kj)*winv(jnl,ki)-winv(jnl,kj)*winv(i,ki)
         tablec=cdet*(tabpipSz(i)/tabpipSz(jn))**2*jhop(i,jn)
        elseif(iconf(i).eq.1.and.iconf(jn).eq.-1) then
         ki=kelup(i)
         kj=keldo(jn)
         il=i+L
         cdet=winv(il,kj)*winv(jn,ki)-winv(jn,kj)*winv(il,ki)
         tablec=cdet*(tabpipSz(jn)/tabpipSz(i))**2*jhop(i,jn)
        endif

       endif
      ENDIF

      table=abs(tablec)

      return
      end

      subroutine uptablehubf(L,irange,ivicn,maxmulti,imulti
     >,table,tabpip,tabpipSz,thop,jhop,winv,iconf
     >,kelup,keldo,ioptPH)
      implicit none
      integer*4 L,i,ii,j,jn,ioptPH,ki,kj,il,jnl,maxmulti,irange,icount
      integer*4 ivicn(L,maxmulti,*),imulti(*)
      integer*2 iconf(L),kelup(L),keldo(L)
      real*8 tabpip(*),tabpipSz(*),thop(L,L),jhop(L,L)
      complex(8) table(L,*),winv(2*L,*)

      if(ioptPH.ge.0) then 

       icount=0
       do j=1,irange
        do ii=1,imulti(j)
         icount=icount+1
         do i=1,L
          table(i,icount)=dcmplx(0.d0,0.d0)
          jn=abs(ivicn(i,ii,j))
          IF(jn.eq.0) THEN
           table(i,icount)=dcmplx(0.d0,0.d0)
          ELSE
           if(iconf(i).eq.iconf(jn)) then
            table(i,icount)=dcmplx(0.d0,0.d0)
           elseif(iconf(i).eq.0.and.iconf(jn).eq.-1) then
            kj=keldo(jn)
            il=i+L
            table(i,icount)=winv(il,kj)
     >*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)*thop(i,jn)
            table(i,icount)=-table(i,icount)
           elseif(iconf(i).eq.-1.and.iconf(jn).eq.2) then
            kj=kelup(jn)
            table(i,icount)=winv(i,kj)
     >*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)*thop(i,jn)
           elseif(iconf(i).eq.2.and.iconf(jn).eq.0) then
            ki=kelup(i)
            kj=keldo(i)
            jnl=jn+L
            table(i,icount)=
     >(winv(jn,ki)*winv(jnl,kj)-winv(jnl,ki)*winv(jn,kj))
     >*(tabpipSz(jn)/tabpipSz(i))**2*jhop(i,jn)
           endif
          ENDIF
         enddo
        enddo
       enddo

      else

       icount=0
       do j=1,irange
        do ii=1,imulti(j)
         icount=icount+1
         do i=1,L
          table(i,icount)=dcmplx(0.d0,0.d0)
          jn=abs(ivicn(i,ii,j))
          IF(jn.eq.0) THEN
           table(i,icount)=dcmplx(0.d0,0.d0)
          ELSE
           if(iconf(i).eq.iconf(jn)) then
            table(i,icount)=dcmplx(0.d0,0.d0)
           elseif(iconf(i).eq.0.and.iconf(jn).eq.1) then
            kj=kelup(jn)
            table(i,icount)=winv(i,kj)
     >*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)*thop(i,jn)
           elseif(iconf(i).eq.0.and.iconf(jn).eq.-1) then
            kj=keldo(jn)
            il=i+L
            table(i,icount)=winv(il,kj)
     >*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)*thop(i,jn)
           elseif(iconf(i).eq.-1.and.iconf(jn).eq.1) then
            kj=kelup(jn)
            ki=keldo(i)
            jnl=jn+L
            table(i,icount)=
     >(winv(i,kj)*winv(jnl,ki)-winv(jnl,kj)*winv(i,ki))
     >*(tabpipSz(i)/tabpipSz(jn))**2*jhop(i,jn)
           endif
          ENDIF
         enddo
        enddo
       enddo

      endif

      return
      end

c====================================================================
c====================================================================
c====================================================================

      subroutine branchingo(nw,wconfn,weight,zeta,icdiff,ipip,jbra)
      implicit none 
      integer*4 nw,i,j,ind,icdiff,ni
      integer*4 ipip(*),jbra(nw)
      real*8 weight,try,tryp,dble,dstep
      real*8 zeta(nw+1),wconfn(*)

      weight=wconfn(1)
      do i=2,nw
       weight=weight+wconfn(i)
      enddo

      dstep=zeta(nw+1)/dble(nw)
      do i=1,nw
       zeta(i)=weight*(dstep+(i-1)/dble(nw))
      enddo 
      zeta(nw+1)=weight+1.d0

      ind=1
      try=0.d0
      icdiff=0
      do i=1,nw
       tryp=try+wconfn(i)
       ni=0    
       dowhile(zeta(ind).lt.tryp.and.zeta(ind).ge.try) 
        jbra(ind)=i
        ind=ind+1
        ni=ni+1
       enddo
       try=try+wconfn(i)
       if(ni.ne.0) icdiff=icdiff+1
      enddo

      call upjbra(nw,jbra,ipip,ipip(nw+1))

      do j=1,nw
       wconfn(j)=1.d0
      enddo

      weight=weight/dble(nw)

      return
      end

      subroutine upjbra(nw,jbra,ipip,indz)
      implicit none
      integer*4 j,i,nw,ind        
      integer*4 jbra(*),ipip(*),indz(*)

      do j=1,nw
       ipip(j)=0
      enddo
      do j=1,nw
       ind=jbra(j)
       ipip(ind)=ipip(ind)+1
      enddo
      ind=0
      do j=1,nw
       if(ipip(j).eq.0) then
        ind=ind+1
        indz(ind)=j
       endif
      enddo
      ind=0
      do j=1,nw
       if(ipip(j).ne.0) then
        jbra(j)=j
        do i=2,ipip(j)
         ind=ind+1
         jbra(indz(ind))=j
        enddo
       endif
      enddo

      return
      end

      subroutine reshuffhub(L,Lz,nw,jbra,iconf
     >,table,tabpip,tabpipSz,winv,nelm2,wsto
     >,diagfn,diag,kelup,keldo,enert)
      implicit none
      integer*4 L,nw,j,k,ind,Lz,nelm2
      integer*4 jbra(*)
      integer*2 iconf(L,*),kelup(L,*),keldo(L,*)
      real*8 tabpip(L,*),tabpipSz(L,*),diagfn(*),wsto(*)
     >,diag(*)
      complex(8) table(Lz,*),winv(nelm2,*),enert(*)

      do j=1,nw
       ind=jbra(j)
       if(j.ne.ind) then
        diag(j)=diag(ind)
        diagfn(j)=diagfn(ind)
        enert(j)=enert(ind)
        wsto(j)=wsto(ind)
        call dcopy(L,tabpip(1,ind),1,tabpip(1,j),1)
        call dcopy(L,tabpipSz(1,ind),1,tabpipSz(1,j),1)
        call zcopy(Lz,table(1,ind),1,table(1,j),1)
        call zcopy(nelm2,winv(1,ind),1,winv(1,j),1)
        do k=1,L
         iconf(k,j)=iconf(k,ind)
         kelup(k,j)=kelup(k,ind)
         keldo(k,j)=keldo(k,ind)
        enddo
       endif
      enddo

      return
      end

c====================================================================

      subroutine reweight0(nw,np,npm,etotw,ipip,psip,alpha,sov
     >,sovdt,econf,epst,iweight,L,enert,scalpar)
      implicit none 
      integer*4 nw,np,npp,i,j,k,info,irank,npm,iweight,ind,L
      integer*4 ipip(*)
      real*8 dble,cost,epst
      real*8 psip(*),alpha(*),sov(npm,npm,*),sovdt(*)
     >,etotw(*),scalpar(*)
      complex(8) econf(nw,*),enert(*)

      npp=np+1

      IF(np.ne.0) THEN 

       cost=1.d0/nw
       do i=1,npp
        do j=i,npp
         sov(i,j,1)=dconjg(econf(1,i))*econf(1,j)*cost
        enddo
       enddo
       do k=2,nw
        do i=1,npp
         do j=i,npp
          sov(i,j,1)=sov(i,j,1)+dconjg(econf(k,i))*econf(k,j)*cost
         enddo
        enddo
       enddo

       do i=1,npp
        do j=i+1,npp
         sov(j,i,1)=sov(i,j,1)
        enddo
       enddo

       ind=0
       do j=1,npp
        do i=1,npp
         ind=ind+1
         sovdt(ind)=sovdt(ind)+sov(i,j,1)
         sov(i,j,2)=sovdt(ind)
        enddo
       enddo

       do j=1,npp
        do i=1,nw
         etotw(j)=etotw(j)-cost*L*enert(i)*dconjg(econf(i,j))
        enddo
       enddo

       if(iweight.eq.1) then 

        cost=1.d0/sov(npp,npp,2)

        do i=1,np
         psip(i)=cost*etotw(i)-cost**2*etotw(npp)*sov(i,npp,2)
        enddo
        do i=1,np
         do j=1,np
          sov(i,j,2)=cost*sov(i,j,2)-cost**2*sov(i,npp,2)*sov(j,npp,2)
          sov(i,j,1)=sov(i,j,2)
         enddo
        enddo

        irank=0
        call dgelssn(sov,np,npm,psip,psip(np+1)
     >,epst,irank,ipip,info,0,scalpar)

        cost=0.d0
        do i=1,np
         do j=1,np
          cost=cost+scalpar(i)*scalpar(j)*psip(i)*psip(j)*sov(j,i,1)
         enddo
        enddo

        cost=dsqrt(cost)
        write(6,*) 'Norm correction=',cost

        psip(npp)=cost

        if(info.ne.0) then 
         write(6,*) 'Info=',info
         stop
        else
         do i=1,npp
          alpha(i)=psip(i)
         enddo
        endif

       endif

      ENDIF

      return
      end

      subroutine dgelssn(sov,nmat,npm,alpha,psip,epsdgel,irank
     >,ipsip,info,ip,scalpar)
      implicit none
      integer*4 ip,nmat,i,j,npm,n1,n2,info,irank,lwork,indi,indj
      integer*4 ipsip(*)
      real*8 epsdgel,epsmach
      real*8 sov(npm,npm,*),alpha(*),psip(*),scalpar(*)

      epsmach=1d-14

      n1=nmat+1
      n2=n1+nmat
      lwork=10*nmat

      IF(ip.eq.0) THEN

       do i=1,nmat
        ipsip(i)=1
       enddo          
       do i=1,nmat
        if(abs(sov(i,i,1)*scalpar(i)).gt.epsmach) then 
         psip(n1+i-1)=1.d0/dsqrt(dabs(sov(i,i,1)))
        else
         psip(n1+i-1)=0.d0
         ipsip(i)=0
        endif
       enddo
       indi=0
       do i=1,nmat
        if(ipsip(i).eq.1) then 
         indi=indi+1
         indj=0
         do j=1,nmat
          if(ipsip(j).eq.1) then 
           indj=indj+1
           sov(indi,indj,2)=sov(i,j,1)*psip(n1+i-1)*psip(n1+j-1)
          endif 
         enddo
        endif
       enddo

       call dsyev('N','L',indi,sov(1,1,2),npm,psip,psip(n2)
     >,lwork,info)

       write(6,*) ' Eigenvalues sov =',(psip(i),i=1,indi)

       if(epsdgel.ne.0.d0) then 
        do i=1,nmat
         sov(i,i,1)=sov(i,i,1)+epsdgel
        enddo
       endif 
 
       irank=nmat

       call mscale(nmat,npm,sov,alpha,scalpar,psip,ipsip)

      ELSE

       irank=0
       call dgetrf(nmat,nmat,sov,npm,ipsip,info)
       call dgetrs('N',nmat,1,sov,npm,ipsip,alpha,nmat,info)
       if(info.eq.0) irank=nmat

      ENDIF

      return
      end

      subroutine mscale(npr,npm,sov,force,scalpar,psip,ipsip)
      implicit none 
      integer*4 np,npm,npr,i,j,k,n1,n2,indi,indj,ind,ii,jj
      integer*4 ipsip(*)
      real*8 cost
      real*8 sov(npm,npm,*),force(*),scalpar(*),psip(*)

c Ordering according to scalpar

      indi=0
      do i=1,npr
       if(ipsip(i).eq.1) then 
        indi=indi+1
        psip(indi)=scalpar(i)
       endif
      enddo

      np=indi

      n1=npr+1
      n2=np+n1

      call dsortx(psip,1,np,ipsip(n1))

c Load Sov(i,j,1) without preconditioning
      ii=0
      do i=1,npr
       if(ipsip(i).eq.1) then 
        ii=ii+1
        jj=0
        do j=1,npr
         if(ipsip(j).eq.1) then 
          jj=jj+1
          sov(ii,jj,2)=sov(i,j,1)
         endif 
        enddo
       endif 
      enddo

      do i=1,np 
       call dcopy(np,sov(1,i,2),1,sov(1,i,1),1)
      enddo

c Ordering Sov(i,j,2) according to scalpar
      do i=1,np
       indi=ipsip(n2-i)
       do j=1,np
        indj=ipsip(n2-j)
        sov(i,j,2)=sov(indi,indj,1)
       enddo
      enddo

c  Begin Graham-Schmidt
      call dscalzero(np,0.d0,sov,1)
      psip(1)=dsqrt(sov(1,1,2))
      cost=1.d0/psip(1)
      sov(1,1,1)=cost
 
      do i=2,np
       call dgemv('T',i-1,i-1,1.d0,sov,npm,sov(1,i,2),1,0.d0
     >,psip(n1),1)

       call dscalzero(np,0.d0,psip(n2),1)
       psip(n2+i-1)=1.d0
       call dgemv('N',i-1,i-1,-1.d0,sov,npm,psip(n1),1,1.d0,psip(n2),1)
       psip(i)=0.d0
       do j=1,i
        do k=1,i
         psip(i)=psip(i)+psip(n2+j-1)*psip(n2+k-1)*sov(j,k,2)
        enddo
       enddo
       psip(i)=dsqrt(psip(i))
       cost=1.d0/psip(i)
       call dscal(i,cost,psip(n2),1)
       call dscalzero(np-i,0.d0,sov(i+1,i,1),1)
       call dcopy(i,psip(n2),1,sov(1,i,1),1)
      enddo
c  End Graham-Schmidt

      indi=0
      do i=1,npr
       if(ipsip(i).eq.1) then 
        indi=indi+1
        psip(indi)=scalpar(i)
        psip(n2+indi-1)=force(i)
       endif
      enddo

      do i=1,np
       ind=ipsip(n2-i)
       psip(n1+i-1)=psip(ind)  ! ordered scalpar
      enddo

      do i=1,np
       ind=ipsip(n2-i)
       psip(i)=psip(n2+ind-1)  ! ordered force
      enddo

      call dgemv('T',np,np,1.d0,sov,npm,psip,1,0.d0,psip(n2),1)

      do i=1,np
       psip(i)=psip(n2+i-1)*psip(n1+i-1)
      enddo

c Go back to the original basis
      do i=1,np
       psip(n2+i-1)=psip(1)*sov(i,1,1)
      enddo

      do i=2,np
       call daxpy(np,psip(i),sov(1,i,1),1,psip(n2),1)
      enddo

      do i=1,np
       ind=ipsip(n2-i)
       psip(ind)=psip(n2+i-1)
      enddo

      indi=0
      do i=1,npr
       if(ipsip(i).eq.1) then 
        indi=indi+1
        force(i)=psip(indi)
       else
        force(i)=0.d0
       endif
      enddo

c Restore matrix sov for the calculation of the norm
      do i=1,np
       indi=ipsip(n2-i)
       do j=1,np
        indj=ipsip(n2-j)
        sov(indi,indj,1)=sov(i,j,2)
       enddo
      enddo

      do i=1,np 
       call dcopy(np,sov(1,i,1),1,sov(1,i,2),1)
      enddo
      ii=0
      do i=1,npr
       if(ipsip(i).eq.1) then
        ii=ii+1
        jj=0
        do j=1,npr
         if(ipsip(j).eq.1) then
          jj=jj+1
          sov(i,j,1)=sov(ii,jj,2)
         else
          sov(i,j,1)=0.d0 
         endif
        enddo
       else
        do j=1,npr
         sov(i,j,1)=0.d0
        enddo
       endif
      enddo

      return
      end

      subroutine dscalzero(n,zero,vet,m)
      implicit none 
      integer*4 n,m,i
      real*8 zero,vet(m,*) 

      do i=1,n
       vet(1,i)=zero
      enddo

      return
      end 

      subroutine makesovop(nh,nel,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,iesBCS,ieskin,iesmu,iesonsite,iesaf,matrix,ioptPH
     >,z,eig,psip,sov)
      implicit none 
      integer*4 nh,nel,k,i,j,ioptPH,maxmultis
      integer*4 ieskin,iesBCS,iesmu,iesonsite,iesaf,ndim,np
      integer*4 ivicns(nh,maxmultis,*),imultis(*)
      real*8 eig(*)
      complex(8) matrix(2*nh,2*nh),psip(2*nh,*),z(2*nh,*)
     >,sov(2*nh,2*nh,*)
     >,phaserKIN(nh,maxmultis,*),phaserBCS(nh,maxmultis,*)
     >,phaseiKIN(nh,maxmultis,*),phaseiBCS(nh,maxmultis,*)
     >,alphai(nh)

      ndim=2*nh
      np=2*iesBCS+2*ieskin+iesmu+2*iesonsite+iesaf

      do k=1,np

       if(k.le.iesBCS) then 
        call upmatrix(nh,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,k,0,0,0,0,0,0,0,matrix,ioptPH)

       elseif(k.le.2*iesBCS) then 
        call upmatrix(nh,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,0,k-iesBCS,0,0,0,0,0,0,matrix,ioptPH)

       elseif(k.le.2*iesBCS+ieskin) then 
        call upmatrix(nh,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,0,0,k-2*iesBCS,0,0,0,0,0,matrix,ioptPH)

       elseif(k.le.2*iesBCS+2*ieskin) then 
        call upmatrix(nh,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,0,0,0,k-2*iesBCS-ieskin,0,0,0,0,matrix,ioptPH)

       elseif(k.le.2*iesBCS+2*ieskin+iesmu) then 
        call upmatrix(nh,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,0,0,0,0,iesmu,0,0,0,matrix,ioptPH)

       elseif(k.le.2*iesBCS+2*ieskin+iesmu+iesonsite) then 
        call upmatrix(nh,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,0,0,0,0,0,iesonsite,0,0,matrix,ioptPH)

       elseif(k.le.2*iesBCS+2*ieskin+iesmu+2*iesonsite) then 
        call upmatrix(nh,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,0,0,0,0,0,0,iesonsite,0,matrix,ioptPH)

       elseif(k.le.2*iesBCS+2*ieskin+iesmu+2*iesonsite+iesaf) then 
        call upmatrix(nh,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,0,0,0,0,0,0,0,iesaf,matrix,ioptPH)

       endif

       call zgemm('N','N',ndim,ndim,ndim,(1.d0,0.d0),matrix,ndim,z
     >,ndim,(0.d0,0.d0),psip,ndim)
       call zgemm('C','N',ndim,ndim,ndim,(1.d0,0.d0),z,ndim,psip
     >,ndim,(0.d0,0.d0),matrix,ndim)

       do i=1,ndim
        do j=1,ndim
         if(i.gt.nel.and.j.le.nel) then 
          matrix(i,j)=-matrix(i,j)/(eig(i)-eig(j))
         else
          matrix(i,j)=0.d0
         endif
        enddo
       enddo

       do i=nel+1,ndim
        do j=1,nel
         matrix(j,i)=dconjg(matrix(i,j))
        enddo
       enddo

       call zgemm('N','C',ndim,ndim,ndim,(1.d0,0.d0),matrix,ndim,z
     >,ndim,(0.d0,0.d0),psip,ndim)
       call zgemm('N','N',ndim,ndim,ndim,(1.d0,0.d0),z,ndim,psip
     >,ndim,(0.d0,0.d0),sov(1,1,k),ndim)

      enddo

      return
      end

      subroutine upmatrix(nh,ivicns,maxmultis,imultis
     >,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS,alphai
     >,iesBCS_r,iesBCS_i,ieskin_r,ieskin_i
     >,iesmu,iesonsite_r,iesonsite_i,iesaf,matrix,ioptPH)
      implicit none 
      integer*4 nh,ip,iq,k,jn
      integer*4 iesBCS_r,iesBCS_i,ieskin_r,ieskin_i,iesmu
     >,iesonsite_r,iesonsite_i,iesaf,ioptPH,maxmultis
      integer*4 ivicns(nh,maxmultis,*),imultis(*)
      complex(8) ai
      complex(8) matrix(2*nh,2*nh)
     >,phaserKIN(nh,maxmultis,*),phaserBCS(nh,maxmultis,*)
     >,phaseiKIN(nh,maxmultis,*),phaseiBCS(nh,maxmultis,*)
     >,alphai(nh)

      ai=dcmplx(0.d0,1.d0)

      do ip=1,2*nh
       do iq=1,2*nh
        matrix(ip,iq)=dcmplx(0.d0,0.d0)
       enddo
      enddo

      if(iesmu.ne.0) then 
       do ip=1,nh
        matrix(ip,ip)=-1.d0
        if(ioptPH.ge.0) then
         matrix(ip+nh,ip+nh)=1.d0
        else
         matrix(ip+nh,ip+nh)=-1.d0
        endif
       enddo
      endif

      if(ieskin_r.ne.0) then 
       do k=1,imultis(ieskin_r)
        do ip=1,nh
         jn=ivicns(ip,k,ieskin_r)
         if(jn.ne.0) then 
          if(jn.gt.0) then
           matrix(ip,jn)=matrix(ip,jn)+phaserKIN(ip,k,ieskin_r)
           if(ioptPH.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >                         +phaserKIN(ip,k,ieskin_r)
           else
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >                         -dconjg(phaserKIN(ip,k,ieskin_r))
           endif
          elseif(jn.lt.0) then
           jn=-jn
           matrix(ip,jn)=matrix(ip,jn)-phaserKIN(ip,k,ieskin_r)
           if(ioptPH.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >                         -phaserKIN(ip,k,ieskin_r)
           else
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >                         +dconjg(phaserKIN(ip,k,ieskin_r))
           endif
          endif
         endif
        enddo
       enddo
      endif
 
      if(ieskin_i.ne.0) then 
       do k=1,imultis(ieskin_i)
        do ip=1,nh
         jn=ivicns(ip,k,ieskin_i)
         if(jn.ne.0) then 
          if(jn.gt.0) then
           matrix(ip,jn)=matrix(ip,jn)+ai*phaseiKIN(ip,k,ieskin_i)
           if(ioptPH.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >                  +ai*phaseiKIN(ip,k,ieskin_i)
           else
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >                  +ai*dconjg(phaseiKIN(ip,k,ieskin_i))
           endif
          elseif(jn.lt.0) then
           jn=-jn
           matrix(ip,jn)=matrix(ip,jn)-ai*phaseiKIN(ip,k,ieskin_i)
           if(ioptPH.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >                  -ai*phaseiKIN(ip,k,ieskin_i)
           else
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)
     >                  -ai*dconjg(phaseiKIN(ip,k,ieskin_i))
           endif
          endif
         endif
        enddo
       enddo
      endif
 
      if(iesonsite_r.ne.0) then 
       do ip=1,nh
        matrix(ip,ip+nh)=-1.d0
        matrix(ip+nh,ip)=-1.d0
       enddo
      endif

      if(iesonsite_i.ne.0) then 
       do ip=1,nh
        matrix(ip,ip+nh)=ai
        matrix(ip+nh,ip)=-ai
       enddo
      endif

      if(iesBCS_r.ne.0) then 
       do k=1,imultis(iesBCS_r)
        do ip=1,nh
         jn=ivicns(ip,k,iesBCS_r)
         if(jn.gt.0) then 
          matrix(ip+nh,jn)=matrix(ip+nh,jn)
     >                    +phaserBCS(ip,k,iesBCS_r)
          matrix(ip,jn+nh)=matrix(ip,jn+nh)
     >                    +dconjg(phaserBCS(ip,k,iesBCS_r))
         elseif(jn.lt.0) then
          jn=-jn
          matrix(ip+nh,jn)=matrix(ip+nh,jn)
     >                    -phaserBCS(ip,k,iesBCS_r)
          matrix(ip,jn+nh)=matrix(ip,jn+nh)
     >                    -dconjg(phaserBCS(ip,k,iesBCS_r))
         endif
        enddo
       enddo
      endif

      if(iesBCS_i.ne.0) then 
       do k=1,imultis(iesBCS_i)
        do ip=1,nh
         jn=ivicns(ip,k,iesBCS_i)
         if(jn.gt.0) then 
          matrix(ip+nh,jn)=matrix(ip+nh,jn)
     >                    +ai*phaseiBCS(ip,k,iesBCS_i)
          matrix(ip,jn+nh)=matrix(ip,jn+nh)
     >                    -ai*dconjg(phaseiBCS(ip,k,iesBCS_i))
         elseif(jn.lt.0) then
          jn=-jn
          matrix(ip+nh,jn)=matrix(ip+nh,jn)
     >                    -ai*phaseiBCS(ip,k,iesBCS_i)
          matrix(ip,jn+nh)=matrix(ip,jn+nh)
     >                    +ai*dconjg(phaseiBCS(ip,k,iesBCS_i))
         endif
        enddo
       enddo
      endif

      if(iesaf.eq.1) then
       do ip=1,nh
        if(ioptPH.lt.0) then
         matrix(ip,ip+nh)=matrix(ip,ip+nh)+alphai(ip)
         matrix(ip+nh,ip)=matrix(ip+nh,ip)+dconjg(alphai(ip))
        endif
       enddo
      endif

      return
      end

      subroutine upallcorr(nw,nh2,np,sov,winv
     >,kelup,keldo,econf)
      implicit none 
      integer*4 nw,nh,nh2,np,i,j,nhs
      integer*2 kelup(*),keldo(*)
      complex(8) sov(nh2,nh2,*),winv(nh2,*),econf(nw,*)

      nh=nh2/2
      nhs=nh2**2

      do i=1,np
       econf(1,i)=dcmplx(0.d0,0.d0)
      enddo

      do j=1,nh
       if(kelup(j).ne.0) then  
        call zgemv('C',nh2,np,(1.d0,0.d0),sov(1,j,1),nhs
     >,winv(1,kelup(j)),1,(1.d0,0.d0),econf,nw)
       endif
       if(keldo(j).ne.0) then 
        call zgemv('C',nh2,np,(1.d0,0.d0),sov(1,j+nh,1),nhs
     >,winv(1,keldo(j)),1,(1.d0,0.d0),econf,nw)
       endif
      enddo

      return
      end

      subroutine upvpot(imis,nw,L,ivicn,maxmulti,imulti,iconf
     >,econf,ioptPH)
      implicit none
      integer*4 L,i,j,k,jn,imis,nw,ioptPH,maxmulti
      integer*4 ivicn(L,maxmulti,*),imulti(*)
      integer*2 iconf(*)
      real*8 dnn
      complex(8) econf(nw,*)

      if(ioptPH.ge.0) then 
       do k=1,imis 
        dnn=0.d0 
        do i=1,imulti(k)
         do j=1,L
          jn=abs(ivicn(j,i,k))
          if(jn.ne.0) then
           if(iconf(j).ne.2.and.iconf(jn).ne.2) 
     >     dnn=dnn+0.5d0*iconf(j)*iconf(jn)  
          endif
         enddo
        enddo
        econf(1,k)=dnn
       enddo
      else
       do k=1,imis 
        dnn=0.d0 
        do i=1,imulti(k)
         do j=1,L
          jn=abs(ivicn(j,i,k))
          if(jn.ne.0) 
     >    dnn=dnn+0.5d0*abs(iconf(j)*iconf(jn))
         enddo
        enddo
        econf(1,k)=dnn
       enddo
      endif

      return
      end

      subroutine upvpotz(imis,nw,L,ivicn,maxmulti,imulti,iconf
     >,econf,ioptPH)
      implicit none
      integer*4 L,i,j,k,jn,imis,nw,ioptPH,maxmulti
      integer*4 ivicn(L,maxmulti,*),imulti(*)
      integer*2 iconf(*)
      real*8 dnn
      complex(8) econf(nw,*)

      if(ioptPH.ge.0) then 
       do k=1,imis 
        dnn=0.d0 
        do i=1,imulti(k)
         do j=1,L
          jn=abs(ivicn(j,i,k))
          if(jn.ne.0) 
     >    dnn=dnn+0.5d0*abs(iconf(j)*iconf(jn))
         enddo
        enddo
        econf(1,k)=dnn
       enddo
      else
       do k=1,imis
        dnn=0.d0
        do i=1,imulti(k)
         do j=1,L
          jn=abs(ivicn(j,i,k))
          if(jn.ne.0) then
           if(iconf(j).ne.2.and.iconf(jn).ne.2)
     >     dnn=dnn+0.5d0*iconf(j)*iconf(jn)
          endif
         enddo
        enddo
        econf(1,k)=dnn
       enddo
      endif

      return
      end

      subroutine uphop(L,irange,ivicn,imulti,maxmulti
     >,v,vsz,tham,jham,thop,thopv,jhop,jhopv,ioptPH)
      implicit none
      integer*4 L,i,j,k,jn,irange,maxmulti,ioptPH
      integer*4 ivicn(L,maxmulti,*),imulti(*)
      real*8 thop(L,L),thopv(L,L),jhop(L,L),jhopv(L,L)
     >,v(L,L),vsz(L,L),tham(*),jham(*)

      do i=1,L
       do j=1,L
        thop(i,j)=0.d0
        jhop(i,j)=0.d0
        thopv(i,j)=0.d0
        jhopv(i,j)=0.d0
       enddo
      enddo

      do k=1,irange
       do j=1,imulti(k)
        do i=1,L
         jn=ivicn(i,j,k)
         if(jn.gt.0) then 
          thop(i,jn)=tham(k)*v(i,i)/v(i,jn)*vSz(i,i)/vSz(i,jn)
          thopv(i,jn)=v(i,i)/v(i,jn)*vSz(i,i)/vSz(i,jn)
          jhop(i,jn)=0.5d0*jham(k)*(vSz(i,i)/vSz(i,jn))**4
          if(ioptPH.ge.0) jhop(i,jn)=-jhop(i,jn)
          jhopv(i,jn)=(vSz(i,i)/vSz(i,jn))**4
         elseif(jn.lt.0) then
          jn=-jn
          thop(i,jn)=-tham(k)*v(i,i)/v(i,jn)*vSz(i,i)/vSz(i,jn)
          thopv(i,jn)=v(i,i)/v(i,jn)*vSz(i,i)/vSz(i,jn)
          jhop(i,jn)=0.5d0*jham(k)*(vSz(i,i)/vSz(i,jn))**4
          if(ioptPH.ge.0) jhop(i,jn)=-jhop(i,jn)
          jhopv(i,jn)=(vSz(i,i)/vSz(i,jn))**4
         endif 
        enddo
       enddo
      enddo

      return
      end

      subroutine ekdk(nh,irange,ivicn,ivicns,maxmulti,maxmultis
     >,imulti,imultis,phaserKIN,phaserBCS,phaseiKIN,phaseiBCS
     >,iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf,iesslater
     >,vjspin,vjdens,dBCS_r,dBCS_i,dek_r,dek_i,dmu,donsite_r,donsite_i
     >,dAF,alphai,psip,rpsip,matrix,vsz,v,z,muf,tvar,eig,ioptPH)
      implicit none 
      integer*4 nh,iq,ip,info,i,k,j,jn,maxmulti,maxmultis
     >,irange,ioptPH
      integer*4 iesspin,iesdens,iesBCS,ieskin,iesmu,iesonsite,iesaf
     >,iesslater
      integer*4 ivicn(nh,maxmulti,*),ivicns(nh,maxmultis,*)
     >,imulti(*),imultis(*)
      real*8 vjdens(*),vjspin(*),dBCS_r(*),dBCS_i(*),dek_r(*),dek_i(*)
     >,dmu(*),donsite_r(*),donsite_i(*),dAF(*),tvar(*)
      real*8 vsz(nh,nh),v(nh,nh),eig(*),rpsip(*),muf,dek1
      complex(8) ai
      complex(8) matrix(2*nh,2*nh),psip(*),z(2*nh,*)
     >,phaserKIN(nh,maxmultis,*),phaserBCS(nh,maxmultis,*)
     >,phaseiKIN(nh,maxmultis,*),phaseiBCS(nh,maxmultis,*)
     >,alphai(nh)
      COMPLEX(8), dimension(:), allocatable :: work

      if(iesdens.ne.0) then
       do i=1,nh
        do j=1,nh
         v(i,j)=1.d0
        enddo
       enddo
      endif
      if(iesspin.ne.0) then
       do i=1,nh
        do j=1,nh
         vsz(i,j)=1.d0
        enddo
       enddo
      endif

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

c Determinant

      ai=dcmplx(0.d0,1.d0)

      do ip=1,2*nh
       do iq=1,2*nh
        matrix(ip,iq)=dcmplx(0.d0,0.d0)
       enddo
      enddo

      if(iesslater.eq.0) return

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

c Matrix diagonalization

      call zcopy(4*nh*nh,matrix,1,psip,1)

      ALLOCATE(work(20*nh))

      call zheev('V','L',2*nh,psip,2*nh,eig,work,20*nh
     >,rpsip,info)

      do i=1,2*nh
       call zcopy(2*nh,psip(2*nh*(i-1)+1),1,z(1,i),1)
      enddo

      DEALLOCATE(work)

      return
      end

      real*8 function upvham(uvv,jperp,jham,L,indt,ivicn,maxmulti,imulti
     >,iconf,ioptPH)
      implicit none
      integer*4 L,indt,i,j,ii,jn,ioptPH,maxmulti
      integer*4 ivicn(L,maxmulti,*),imulti(*)
      integer*2 iconf(*)
      real*8 uvv,jperp,jham(*)

      if(ioptPH.ge.0) then
       upvham=0.d0
       do i=1,indt
        do ii=1,imulti(i)
         do j=1,L
          jn=abs(ivicn(j,ii,i))
          if(jn.ne.0) then
           if(iconf(j).ne.-1.and.iconf(jn).ne.-1) 
     >     upvham=upvham+uvv*jham(i)*jperp
          endif
         enddo
        enddo
       enddo
      else
       upvham=0.d0
       do i=1,indt
        do ii=1,imulti(i)
         do j=1,L
          jn=abs(ivicn(j,ii,i))
          if(jn.ne.0) then
           if(iconf(j).ne.0.and.iconf(jn).ne.0)
     >     upvham=upvham+uvv*jham(i)*jperp
          endif
         enddo
        enddo
       enddo
      endif

      return
      end

      subroutine upinvhop2(L2,nel,knew,kelk,lnew,kell,winv,psip,epst)
      implicit none
      integer*2 kelk,kell
      integer*4 i,nel,knew,lnew,L2
      real*8 epst
      complex(8) winv(L2,*),c(2,2),psip(L2,*),g

      c(1,1)=winv(lnew,kell)
      c(2,2)=winv(knew,kelk)
      c(1,2)=-winv(knew,kell)
      c(2,1)=-winv(lnew,kelk)

      g=c(1,1)*c(2,2)-c(1,2)*c(2,1)

      if(abs(g).lt.epst) then
       return
      endif

      call zcopy(L2,winv(1,kelk),1,psip,1)
      call zcopy(L2,winv(1,kell),1,psip(1,2),1)
      call zcopy(nel,winv(knew,1),L2,psip(1,3),1)
      psip(kelk,3)=psip(kelk,3)-1.d0
      call zcopy(nel,winv(lnew,1),L2,psip(1,4),1)
      psip(kell,4)=psip(kell,4)-1.d0

      g=-(1.d0,0.d0)/g

      do i=1,nel
       psip(i+nel,3)=c(1,1)*psip(i,3)+c(1,2)*psip(i,4)
       psip(i+nel,4)=c(2,1)*psip(i,3)+c(2,2)*psip(i,4)
      enddo

      call zgemm('N','T',L2,nel,2,g,psip,L2
     >,psip(nel+1,3),L2,(1.d0,0.d0),winv,l2)

      return
      end

C     ==================================================================
      SUBROUTINE DSORTX(COUNT,INUTILE,N,INDEX)
C     ==--------------------------------------------------------------==
C     == Sorting routine for the reciprocal space vectors (g)         ==
C     == KB07AD HANDLES DOUBLE PRECISION VARIABLES                    ==
C     == STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)            ==
C     == THE WORK-SPACE 'MARK' OF LENGTH 50 PERMITS UP TO 2**(50/2)   ==
C     == NUMBERS TO BE SORTED. THIS IS MORE THAN THE IBM VIRTUAL      ==
C     == MEMORY SPACE WILL HOLD .                                     ==
C     ==--------------------------------------------------------------==
      REAL*8 COUNT(N),AV,X
      INTEGER INDEX(N)
      INTEGER INUTILE
      DIMENSION MARK(50)
C     ==--------------------------------------------------------------==
C     ==  SET INDEX ARRAY TO ORIGINAL ORDER .                         ==
C     ==--------------------------------------------------------------==
      DO I=1,N
        INDEX(I)=I
      ENDDO
C     ==--------------------------------------------------------------==
C     == CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED.              ==
C     ==--------------------------------------------------------------==
      IF(N.EQ.1)GOTO 200
      IF(N.GE.1)GOTO 30
      WRITE(6,20)
   20 FORMAT(///20X,65H    KB07AD   NO NUMBERS TO BE SORTED ** RETURN TO
     2 CALLING PROGRAM )
      GOTO 200
C     ==--------------------------------------------------------------==
C     == 'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER  ==
C     == THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.         ==
C     ==--------------------------------------------------------------==
   30 M=12
C     ==--------------------------------------------------------------==
C     == SET UP INITIAL VALUES.                                       ==
C     ==--------------------------------------------------------------==
      LA=2
      IS=1
      IF=N
      DO 190 MLOOP=1,N
C     ==--------------------------------------------------------------==
C     ==  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE. ==
C     ==--------------------------------------------------------------==
        IFKA=IF-IS
        IF((IFKA+1).GT.M)GOTO 70
C     ==--------------------------------------------------------------==
C     == FINAL SORTING  ( A SIMPLE BUBBLE SORT )                      ==
C     ==--------------------------------------------------------------==
        IS1=IS+1
        DO 60 J=IS1,IF
          I=J
   40     IF(COUNT(I-1).LT.COUNT(I))GOTO 60
          IF(COUNT(I-1).GT.COUNT(I))GOTO 50
          IF(INDEX(I-1).LT.INDEX(I))GOTO 60
   50     AV=COUNT(I-1)
          COUNT(I-1)=COUNT(I)
          COUNT(I)=AV
          INT=INDEX(I-1)
          INDEX(I-1)=INDEX(I)
          INDEX(I)=INT
          I=I-1
          IF(I.GT.IS)GOTO  40
   60   CONTINUE
        LA=LA-2
        GOTO 170
C     ==--------------------------------------------------------------==
C     ==                      *  QUICKSORT        **                  ==
C     == SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS  ==
C     == THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S==
C     == HIGHEST ADDRESS.                                             ==
C     ==--------------------------------------------------------------==
   70   IY=(IS+IF)/2
        X=COUNT(IY)
        INTEST=INDEX(IY)
        COUNT(IY)=COUNT(IF)
        INDEX(IY)=INDEX(IF)
C     ==--------------------------------------------------------------==
C     == THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END ==
C     == OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE   ==
C     == OF X .                                                       ==
C     ==--------------------------------------------------------------==
        K=1
        IFK=IF
C     ==--------------------------------------------------------------==
C     == WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE ==
C     == INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS   ==
C     == NECESSARY, UNTIL THEY MEET .                                 ==
C     ==--------------------------------------------------------------==
        DO 110 I=IS,IF
          IF(X.GT.COUNT(I))GOTO 110
          IF(X.LT.COUNT(I))GOTO 80
          IF(INTEST.GT.INDEX(I))GOTO 110
   80     IF(I.GE.IFK)GOTO 120
          COUNT(IFK)=COUNT(I)
          INDEX(IFK)=INDEX(I)
          K1=K
          DO 100 K=K1,IFKA
            IFK=IF-K
            IF(COUNT(IFK).GT.X)GOTO 100
            IF(COUNT(IFK).LT.X)GOTO 90
            IF(INTEST.LE.INDEX(IFK))GOTO 100
   90       IF(I.GE.IFK)GOTO 130
            COUNT(I)=COUNT(IFK)
            INDEX(I)=INDEX(IFK)
            GO TO 110
  100     CONTINUE
          GOTO 120
  110   CONTINUE
C     ==--------------------------------------------------------------==
C     == RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER  ==
C     == WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO ==
C     == 2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL ==
C     == TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED==
C     == INDEPENDENTLY .                                              ==
C     ==--------------------------------------------------------------==
  120   COUNT(IFK)=X
        INDEX(IFK)=INTEST
        IP=IFK
        GOTO 140
  130   COUNT(I)=X
        INDEX(I)=INTEST
        IP=I
C     ==--------------------------------------------------------------==
C     ==  STORE THE LONGER SUBDIVISION IN WORKSPACE.                  ==
C     ==--------------------------------------------------------------==
  140   IF((IP-IS).GT.(IF-IP))GOTO 150
        MARK(LA)=IF
        MARK(LA-1)=IP+1
        IF=IP-1
        GOTO 160
  150   MARK(LA)=IP-1
        MARK(LA-1)=IS
        IS=IP+1
C     ==--------------------------------------------------------------==
C     == FIND THE LENGTH OF THE SHORTER SUBDIVISION.                  ==
C     ==--------------------------------------------------------------==
  160   LNGTH=IF-IS
        IF(LNGTH.LE.0)GOTO 180
C     ==--------------------------------------------------------------==
C     == IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE==
C     ==--------------------------------------------------------------==
        LA=LA+2
        GOTO 190
  170   IF(LA.LE.0)GOTO 200
C     ==--------------------------------------------------------------==
C     == OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT==
C     ==--------------------------------------------------------------==
  180   IF=MARK(LA)
        IS=MARK(LA-1)
  190 CONTINUE
C     ==--------------------------------------------------------------==
  200 RETURN
      END
C     ==================================================================
