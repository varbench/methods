      program turbohub      
      implicit none 
      integer*4 nx,ny,idim,ivc,indt,izeth,L,nelm,nwm,npm,nelm2
     >,iscramax,L2,Lz,nw,nel,ind1,ind2,ngr,ngsr,nelupr,neldor
     >,nelup,neldo,inext,ibinit,itest,itestr,ioptir
     >,ngg,ngn,ngs,ng,ngen,nscra,nbra,nweight,iseed,iopt,iread
     >,iesm,iesd,iesh,iback(4),ibacktot,iesfree,iessw,ieskin,iesup
     >,iesslater,np,nmat,npsov,nmats,nwm2,nwdim,nwhd,nwfree
     >,nwsw,nwkin,nwup,indc,icdiff,iend,iendr,nwr,Lr,nbrar,npr,nwback
     >,ireadr,nacc,mclock,nleft,ntry,i,j,k,ii,kk,iout,jout,indvic
     >,indspin,info,nint,ioutup,joutup,ioutdo,joutdo,ioutr,cont,nn
     >,nraw,iflag,ind,stripe
      parameter(nx=48,L=288,idim=2,ivc=2,izeth=ivc*idim,indt=2*izeth
     >,nelm=L,nwm=2,npm=180,npsov=npm*npm)
      parameter(iscramax=8*L*L+22*L+npsov+13*npm) 
      real*8 etry,tbra,gamma,epst,epsdgel,tau,U
     >,lambda,muf,etryr,tbrar,wbra,drand1
     >,sumdiff,timescra,timep,veff,tleft,costw,pdiag,ttry
     >,ratio,ener,ener2,diagsum,wtot,rata,time
     >,sflip,p1,p2,costz,tacc,pi,uniform_charge,uniform_spin
     <,uniform_charge_im,uniform_spin_im
      integer*2 ivic(L,indt),ivicr(L,indt),icord(L,idim),iconf(L,nwm)
     >,jbraw(nwm),kelup(L,nwm),keldo(L,nwm),irif(L),multirif(L)
     >,ivicns(L,8,L),ivicn(L,8,L),ivicnr(L,8,L),ivicnb(L,8,L)
     >,kel(L),iconfnew(L)
      integer*4 indtn(L),indtns(L),indtnr(L),indtnb(L),jbra(nwm)
     >,iconfm((nwm*L)/16+1),ipsip(L+2*nwm+2*npm),ipos(npm)
      real*4 wconfw(nwm)
      real*8 v(L,L),vSz(L,L),vhd(L,L),t(izeth),alphav(izeth)
     >,hopupdo(2),z(2*L,2*L),vjz(L),vj(L),vjh(L),ddw(L),dsw(L)
     >,dek(L),dup(L),par(L),matrix(2*L,2*L),eig(2*L)
     >,opmat(2*L,2*L,npm),back(L,4),alpha(npm)
     >,gh(izeth),ghSz(izeth),ghhd(izeth),thop(L,izeth),thopv(L,izeth)
     >,jhop(L,izeth),jhopv(L,izeth),sfhopv(L,izeth)
     >,wconfn(nwm),econf(Nwm*npm)
     >,zeta(nwm+1),tabpip(L,nwm),tabpipSz(L,nwm),tabpiphd(2*L,nwm)
     >,table(L,indt,nwm),diag(nwm),diagfn(nwm),wsto(nwm)
     >,enert(nwm),enert2(nwm),diagpart(nwm)
     >,sovdt(npsov),sov(npm,npm,2),etotw(npm)
     >,zet(nelm,L),cdet(L,L),matr(L,L)
     >,inver(nelm,nelm,nwm),amat(L,L,2),scalpar(npm)
     >,psip(iscramax),work(L),psim(2*L,2*L),zetscra(nelm)
     >,total_charge(0:nx/2),total_spin(0:nx/2)
     >,charge_parity(nx/2),spin_parity(nx/2)
     >,charge_parity_im(nx/2),spin_parity_im(nx/2)
      logical iocc(L)
      character*11 nomefile

      namelist /montecarlo/ itestr,iseed,iopt,iread,nw,ngen,nscra,epst
      namelist /fixednode/ etry,tbra,gamma
      namelist /vmc/ nbra,sflip
      namelist /optimization/ nbra,sflip,nweight,ibinit,epsdgel,tau,nraw

      nomefile='randomseeds' !//char(0)

      L2=2*L
      nelm2=nelm*nelm
      Lz=L*indt

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
      sflip=0.d0
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
       write(6,*) 'sflip    =',sflip
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
       write(6,*) 'sflip    =',sflip
       write(6,*) 'nweight  =',nweight
       write(6,*) 'ibinit   =',ibinit
       write(6,*) 'epsdgel  =',epsdgel
       write(6,*) 'tau      =',tau
       write(6,*) 'nraw     =',nraw
      else
       write(6,*) 'Wrong itestr',itestr
       stop
      endif

      do i=1,npm
       scalpar(i)=1.d0
      enddo
      do i=1,nraw
       read(5,*) nn,tacc,(ipos(j),j=1,nn)
       do j=1,nn
        scalpar(ipos(j))=tacc
       enddo
      enddo

      itest=abs(itestr)
      p1=0.5d0*(1.d0-sflip)
      p2=1.d0-sflip

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
         
      read(10) v,vSz,vhd,t,U,ivicr,icord,irif,nelup,neldo
     >,z,ivicn,indtn,ivicns,indtns,ivicnr,indtnr,ivicnb,indtnb
     >,vjz,vj,vjh,back,ddw,dsw,dek,dup,iesm,iesd,iesh,iback,iesfree
     >,iessw,ieskin,iesup,muf,stripe,ioptir,alphav,par,eig

      do i=1,L
       multirif(i)=0
       if(indtnr(i).ne.0) then
        multirif(i)=1
        ind=abs(ivicnr(1,1,i))
        if(irif(ind).ne.ind) multirif(i)=2
       endif
      enddo

      hopupdo(1)=1.d0
      hopupdo(2)=1.d0
      if(ioptir.ge.0) then
       hopupdo(2)=-1.d0
      endif

      do i=1,indt
       do j=1,L
        ivic(j,i)=abs(ivicr(j,i))
       enddo
      enddo

      do i=1,npm
       alpha(i)=0.d0
      enddo

      ibacktot=iback(1)+iback(2)+iback(3)+iback(4)
      np=iesm+iesd+iesh+ibacktot+iesfree+iessw+iesup+ieskin
      nmat=np+1
      nmats=iesfree+iessw+iesup+ieskin
      if(itestr.eq.-2) then
       if(nmat.gt.npm) then
        write(6,*) 'Increase npm',nmat,npm
        stop
       endif
      endif

      if(nraw.ne.0) then
       write(6,*)
       write(6,*) 'Scalpar'
       do i=1,np
        write(6,*) i,scalpar(i)
       enddo
       write(6,*)
      endif

      if(nw.gt.nwm) then
       write(6,*) ' Too many walkers, recompile with nwm>=',nw
       stop
      endif

      nwm2=0
      nwdim=iesm*nw
      nwhd=(iesm+iesd)*nw
      nwback=(iesm+iesd+iesh)*nw
      nwfree=(iesm+iesd+iesh+ibacktot)*nw
      nwsw=(iesm+iesd+iesh+ibacktot+iesfree)*nw
      nwup=(iesm+iesd+iesh+ibacktot+iesfree+iessw)*nw
      nwkin=(iesm+iesd+iesh+ibacktot+iesfree+iessw+iesup)*nw 

      nel=nelup+neldo

      if(iesfree+iessw+ieskin+iesup.ne.0) then 
       iesslater=1
      else
       iesslater=0
      endif

      if(itestr.eq.-2.and.iesslater.eq.1)
     > call makesovop(L,nelup,neldo,icord,indtn,ivicn,indtns
     >,ivicns,indtnr,ivicnr,iesfree,iessw,iesup,ieskin,matrix,par
     >,stripe,ioptir,z,eig,psip,opmat)

      call uphop(L,izeth,ivicr,v,vsz,vhd
     >,gh,ghsz,ghhd,t,thop,thopv,jhop,jhopv,sfhopv,ioptir)

      if(iopt.eq.1) then

       call initspin(nw,L,nelup,neldo,nel,iconf,kelup,keldo
     >,wconfn,z,ipsip,psip,par,epst,ioptir)

       ng=0
       ngs=0
       iend=0

      elseif(iopt.eq.0.or.iopt.eq.2) then

       call read_seeds(nomefile)
       rewind(11)
       read(11) ngr,ngsr,etryr,nwr,Lr,nbrar,tbrar,npr,ireadr,iendr
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

       if(Lr.ne.L) then
        write(6,*) 'The calculation must continue with L=',Lr
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

       read(11) ((iconf(i,j),i=1,L),j=1,nw),(wconfn(j),j=1,nw)

       do i=1,nw
        wconfn(i)=1.d0
        call upkel(L,nelup,kelup(1,i),keldo(1,i),iconf(1,i))            
       enddo

       rewind(12) 
       if(itestr.eq.-2) rewind(13)
       if(iread.eq.1) rewind(15)
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

      do j=1,Nwm*npm
       econf(j)=1.d0
      enddo
      do j=1,nw
       jbraw(j)=j
       jbra(j)=j
      enddo

      do j=1,npm
       etotw(j)=0.d0
      enddo
      do j=1,npsov
       sovdt(j)=0.d0 
      enddo

      icdiff=0
      wbra=1.d0

      nacc=0
      sumdiff=0.d0

      time=mclock()/100.d0
      timescra=0.d0 

      lambda=-L*etry          

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

        call upscratchub(L,Nw,nelm,nelup,neldo,nel
     >,tabpip,tabpipSz,tabpiphd,diag,v,vSz,vhd,U,t,iconf
     >,kelup,keldo,psip,ipsip,inver,z,ioptir,ivicnb,indtnb
     >,back,iback,zetscra,hopupdo)
       
        if(itest.ne.2.and.i.eq.iend+1) then 
         do j=1,nw
          call uptablehubf(L,idim,ivc,indt,nelm,nel,ivic
     >,table(1,1,j),tabpip(1,j),tabpipSz(1,j),tabpiphd(1,j),thop,jhop
     >,iconf(1,j),iconfnew,z,hopupdo,inver(1,1,j),ipsip,zet
     >,kelup(1,j),keldo(1,j),kel,cdet,matr,ioptir,iocc
     >,ivicnb,indtnb,iback,back)

          call energy(Lz,table(1,1,j),diag(j),wsto(j),veff,enert(j))
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
          iout=drand1()*L+1
          indvic=izeth*drand1()+1
          if(nx.eq.L) indvic=1
          costz=drand1()
          if(costz.lt.p1) then
           indspin=1
          elseif(costz.lt.p2) then
           indspin=-1
          else
           indspin=2 ! Spin flip
          endif
          call ratiovmc(iout,indvic,indspin,L,nelm,nel,ivic
     >,ratio,tabpip(1,j),tabpipSz(1,j),tabpiphd(1,j),thopv,jhopv
     >,z,hopupdo,zet,iconf(1,j),iconfnew,kelup(1,j),keldo(1,j)
     >,ioptir,indt,iocc,sfhopv,inver(1,1,j),ipsip,cont,cdet,matr
     >,kel,ivicnb,indtnb,iback,back)
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

           call random(L,Lz,zeta,table(1,1,j),iout,indvic,gamma)

           if(indvic.le.ivc*idim) then
            indspin=1
           elseif(indvic.gt.ivc*idim) then
            indspin=-1
            indvic=indvic-ivc*idim
           endif
c Here we call ratiovmc to have cdet, zet, and cont
          call ratiofn(iout,indvic,indspin,L,nelm,nel,ivic
     >,z,hopupdo,zet,iconf(1,j),iconfnew,kelup(1,j),keldo(1,j)
     >,ioptir,indt,iocc,inver(1,1,j),ipsip,cont,cdet,matr,kel
     >,ivicnb,indtnb,iback,back)
          endif

          if(indspin.ne.2) then

           if(iconf(iout,j).eq.0.or.
     1        iconf(ivic(iout,indvic),j).eq.0) then
            if(iconf(iout,j).eq.0) then
             jout=ivic(iout,indvic)
            else
             jout=iout
             iout=ivic(jout,indvic)
            endif
           elseif(iconf(iout,j).ne.0.and.
     1            iconf(ivic(iout,indvic),j).ne.0) then
            if(indspin.eq.1) then
             if(iconf(iout,j).eq.-1) then
              jout=ivic(iout,indvic)
             else
              jout=iout
              iout=ivic(jout,indvic)
             endif
            elseif(indspin.eq.-1) then
             if(iconf(iout,j).eq.1) then
              jout=ivic(iout,indvic)
             else
              jout=iout
              iout=ivic(jout,indvic)
             endif
            endif
           endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           if(cont.eq.1)then
            call upinvhop(L,nelm,nel,zet,inver(1,1,j),psip
     >,epst,cdet,kel) 
           else
            call upinvhopbig(L,nelm,nel,zet,inver(1,1,j)
     >,epst,cdet,kel,cont,psip,work,ipsip)
           endif        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           if(indspin.eq.1)then
       call uptabtothubf(L,nel,tabpip(1,j),tabpipSz(1,j),tabpiphd(1,j)
     >,iout,jout,iconf(1,j),v(1,iout),v(1,jout),vSz(1,iout),vSz(1,jout)
     >,vhd(1,iout),vhd(1,jout),kelup(1,j),1,ioptir)
           else
       call uptabtothubf(L,nel,tabpip(1,j),tabpipSz(1,j),tabpiphd(1,j)
     >,iout,jout,iconf(1,j),v(1,iout),v(1,jout),vSz(1,iout),vSz(1,jout)
     >,vhd(1,iout),vhd(1,jout),keldo(1,j),-1,ioptir)
           endif


          elseif(indspin.eq.2) then
           if(iconf(iout,j).lt.iconf(ivic(iout,indvic),j)) then
            jout=ivic(iout,indvic)
           else
            jout=iout
            iout=ivic(jout,indvic)
           endif
           if(ioptir.ge.0) then
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

           ioutr=ioutdo+L

       if(cont.eq.2) then 
        call upinvhop2(L,nelm,nel,zet,inver(1,1,j),psip,ipsip
     >,epst,cdet,kel)
       elseif(cont.gt.2) then
        call upinvhopbig(L,nelm,nel,zet,inver(1,1,j)
     >,epst,cdet,kel,cont,psip,work,ipsip)
       endif

       call uptabtothubf(L,nel,tabpip(1,j),tabpipSz(1,j),tabpiphd(1,j)
     >,ioutup,joutup,iconf(1,j),v(1,ioutup),v(1,joutup),vSz(1,ioutup)
     >,vSz(1,joutup),vhd(1,ioutup),vhd(1,ioutup),kelup(1,j),1,ioptir)
       call uptabtothubf(L,nel,tabpip(1,j),tabpipSz(1,j),tabpiphd(1,j)
     >,ioutdo,joutdo,iconf(1,j),v(1,ioutdo),v(1,joutdo),vSz(1,ioutdo)
     >,vSz(1,joutdo),vhd(1,ioutdo),vhd(1,joutdo),keldo(1,j),-1,ioptir)

          endif

          diag(j)=0.d0
          if(ioptir.ge.0) then 
           do ii=1,L
            if(iconf(ii,j).eq.1)  diag(j)=diag(j)-u 
           enddo
          else
           do ii=1,L
            if(iconf(ii,j).eq.2)  diag(j)=diag(j)-u 
           enddo
          endif

          if(itest.eq.1) then 
           call uptablehubf(L,idim,ivc,indt,nelm,nel,ivic
     >,table(1,1,j),tabpip(1,j),tabpipSz(1,j),tabpiphd(1,j),thop,jhop
     >,iconf(1,j),iconfnew,z,hopupdo,inver(1,1,j),ipsip,zet
     >,kelup(1,j),keldo(1,j),kel,cdet,matr,ioptir,iocc
     >,ivicnb,indtnb,iback,back)

           call energy(Lz,table(1,1,j),diag(j),wsto(j),veff,enert(j))
           diagfn(j)=diag(j)+(1.d0+gamma)*veff
          endif

         endif    ! fine if accepted meas

        enddo ! end do for the tleft or nleft

       enddo ! end do for the walkers

c Calculation of observables
       do j=1,nw
        call uptablehubf(L,idim,ivc,indt,nelm,nel,ivic
     >,table(1,1,j),tabpip(1,j),tabpipSz(1,j),tabpiphd(1,j),thop,jhop
     >,iconf(1,j),iconfnew,z,hopupdo,inver(1,1,j),ipsip,zet
     >,kelup(1,j),keldo(1,j),kel,cdet,matr,ioptir,iocc
     >,ivicnb,indtnb,iback,back)

        IF(itestr.eq.-2) THEN
         if(iesm.ne.0) then
          call upvpotz(iesm,nw,L,indtnr,ivicnr,multirif,iconf(1,j)
     >,econf(nwm2+j),ioptir)
         endif
         if(iesd.ne.0) then 
          call upvpot(iesd,nw,L,indtnr,ivicnr,multirif,iconf(1,j)
     >,econf(nwdim+j),ioptir)
         endif
         if(iesh.ne.0) then 
          call upvpothd(iesh,nw,L,indtnr,ivicnr,multirif,iconf(1,j)
     >,econf(nwhd+j),ioptir)
         endif

         if(ibacktot.ne.0) then
          call upback(iback,nw,L,L2,nelm,nel
     >,iconf(1,j),indtnb,ivicnb,kelup(1,j),keldo(1,j),inver(1,1,j)
     >,z,hopupdo,psip,econf(nwback+j),ioptir,iocc)
         endif
         if(iesslater.eq.1) then
          call upallcorr(L,nw,L2,nel,nelm,nelup,nmats
     >,opmat,inver,z,amat,indtnb,ivicnb,iconf(1,j),hopupdo,back,psim
     >,kelup(1,j),keldo(1,j),econf(j+nwfree),ioptir,iback)
         endif
         

         econf(np*nw+j)=1.d0

        ENDIF

        call energy(Lz,table(1,1,j),diag(j),wsto(j),veff,enert(j))
        enert(j)=enert(j)/dble(L)

        enert2(j)=enert(j)**2
        diagpart(j)=-diag(j)/dble(L)


cccccccc calcolo parita
       ny=L/nx
       pi=2.d0*dasin(1.d0)

         do ii=0,nx/2
          total_charge(ii)=0.d0
          total_spin(ii)=0.d0
         enddo

       if(ioptir.ge.0) then

         do ii=1,nx/2
          total_charge(ii)=total_charge(ii-1)
          total_spin(ii)=total_spin(ii-1)
          do kk=0,ny-1
           if(mod(iconf(ii+kk*nx,j),2).eq.0) then
            total_charge(ii)=total_charge(ii)+1.d0
           else
            total_charge(ii)=total_charge(ii)+iconf(ii+kk*nx,j)+1.d0
           endif
           if(abs(iconf(ii+kk*nx,j)).ne.1) then
            total_spin(ii)=total_spin(ii)+iconf(ii+kk*nx,j)-1.d0
           endif
          enddo
          total_charge(ii)=total_charge(ii)-dble(ny)
     <*dble(nelup+L-neldo)/dble(L)

          charge_parity(ii)=dcos(total_charge(ii)*pi/dble(ny))
          spin_parity(ii)=dcos(total_spin(ii)*pi/dble(ny))
          charge_parity_im(ii)=dsin(total_charge(ii)*pi/dble(ny))
          spin_parity_im(ii)=dsin(total_spin(ii)*pi/dble(ny))

         enddo


        else

         do ii=1,nx/2
          total_charge(ii)=total_charge(ii-1)
          total_spin(ii)=total_spin(ii-1)
          do kk=0,ny-1
           total_charge(ii)=total_charge(ii)+abs(iconf(ii+kk*nx,j))
           if(iconf(ii+kk*nx,j).ne.2) then
            total_spin(ii)=total_spin(ii)+iconf(ii+kk*nx,j)
           endif
          enddo
          total_charge(ii)=total_charge(ii)-dble(ny)
     <*dble(nelup+neldo)/dble(L)

          charge_parity(ii)=dcos(total_charge(ii)*pi/dble(ny))
          spin_parity(ii)=dcos(total_spin(ii)*pi/dble(ny))
          charge_parity_im(ii)=dsin(total_charge(ii)*pi/dble(ny))
          spin_parity_im(ii)=dsin(total_spin(ii)*pi/dble(ny))

         enddo

        endif

        uniform_charge=(charge_parity(nx/2)+charge_parity(nx/2-1))/2.d0
        uniform_spin=(spin_parity(nx/2)+spin_parity(nx/2-1))/2.d0
        uniform_charge_im=(charge_parity_im(nx/2)
     <+charge_parity_im(nx/2-1))/2.d0
        uniform_spin_im=(spin_parity_im(nx/2)
     <+spin_parity_im(nx/2-1))/2.d0



ccccccccccccccccccccccccccccccccccccccc

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
        call branchingo(Nw,wconfn,wbra,zeta,icdiff,ipsip,jbra)
        sumdiff=sumdiff+icdiff

        do j=1,nw
         jbraw(j)=jbra(j)
        enddo

       elseif(itestr.eq.-2) then 

        IF(i.eq.inext) THEN 

         call reweight0(Nw,np,npm,nwm,etotw,ipsip,psip,alpha,sov
     >,sovdt,econf,epsdgel,1,L,enert,scalpar)

         indc=0
         if(iesm.ne.0) then 
          do k=indc+1,indc+iesm
           vjz(k-indc)=vjz(k-indc)+alpha(k)*tau
           alpha(k)=vjz(k-indc)
          enddo
         endif
         indc=indc+iesm
         if(iesd.ne.0) then 
          do k=indc+1,indc+iesd
           vj(k-indc)=vj(k-indc)+alpha(k)*tau
           alpha(k)=vj(k-indc)
          enddo
         endif
         indc=indc+iesd
         if(iesh.ne.0) then 
          do k=indc+1,indc+iesh
           vjh(k-indc)=vjh(k-indc)+alpha(k)*tau
           alpha(k)=vjh(k-indc)
          enddo
         endif
         indc=indc+iesh
         if(iback(1).ne.0)then
          do k=indc+1,indc+iback(1)
           back(k-indc,1)=back(k-indc,1)+alpha(k)*tau
           alpha(k)=back(k-indc,1)
          enddo
         endif
         indc=indc+iback(1)
         if(iback(2).ne.0)then
          do k=indc+1,indc+iback(2)
           back(k-indc,2)=back(k-indc,2)+alpha(k)*tau
           alpha(k)=back(k-indc,2)
          enddo
         endif
         indc=indc+iback(2)
         if(iback(3).ne.0)then
          do k=indc+1,indc+iback(3)
           back(k-indc,3)=back(k-indc,3)+alpha(k)*tau
           alpha(k)=back(k-indc,3)
          enddo
         endif
         indc=indc+iback(3)
         if(iback(4).ne.0)then
          do k=indc+1,indc+iback(4)
           back(k-indc,4)=back(k-indc,4)+alpha(k)*tau
           alpha(k)=back(k-indc,4)
          enddo
         endif
         indc=indc+iback(4)
         if(iesfree.ne.0) then 
          do k=indc+1,indc+iesfree
           ddw(k-indc)=ddw(k-indc)+alpha(k)*tau
           alpha(k)=ddw(k-indc)
          enddo
         endif
         indc=indc+iesfree
         if(iessw.ne.0) then 
          do k=indc+1,indc+iessw
           dsw(k-indc)=dsw(k-indc)+alpha(k)*tau
           alpha(k)=dsw(k-indc)
          enddo
         endif
         indc=indc+iessw
         if(iesup.ne.0) then 
          do k=indc+1,indc+iesup
           dup(k-indc)=dup(k-indc)+alpha(k)*tau
           alpha(k)=dup(k-indc)
          enddo
         endif
         indc=indc+iesup
         if(ieskin.ne.0) then 
          do k=indc+1,indc+ieskin
           dek(k-indc)=dek(k-indc)+alpha(k)*tau
           alpha(k)=dek(k-indc)
          enddo
         endif
         indc=indc+ieskin

         call ekdk(L,nelup,neldo,icord,indtn,ivicn,indtns,ivicns
     >,indtnr,ivicnr,iesm,iesd,iesh,iesfree,iessw,iesup
     >,ieskin,iesslater,vjz,vj,vjh,ddw,dsw,dek,dup,psip,matrix,par
     >,vsz,v,vhd,z,muf,stripe,ivicr,indt,alphav,ioptir,eig)

         if(iesslater.eq.1)
     > call makesovop(L,nelup,neldo,icord,indtn,ivicn,indtns
     >,ivicns,indtnr,ivicnr,iesfree,iessw,iesup,ieskin,matrix,par
     >,stripe,ioptir,z,eig,psip,opmat)

         call uphop(L,izeth,ivicr,v,vsz,vhd
     >,gh,ghsz,ghhd,t,thop,thopv,jhop,jhopv,sfhopv,ioptir)

         call upscratchub(L,Nw,nelm,nelup,neldo,nel
     >,tabpip,tabpipSz,tabpiphd,diag,v,vSz,vhd,U,t,iconf
     >,kelup,keldo,psip,ipsip,inver,z,ioptir,ivicnb,indtnb
     >,back,iback,zetscra,hopupdo)

         write(6,*) ' par   =',i,(alpha(k),k=1,nmat)

         ngn=ngn+1
         write(13) (alpha(k),k=1,nmat)

         do k=1,npm
          etotw(k)=0.d0
         enddo
         do k=1,npsov
          sovdt(k)=0.d0
         enddo

         inext=inext+nweight

        ELSE  ! i.eq.init

         call reweight0(Nw,np,npm,nwm,etotw,ipsip,psip,alpha,sov
     >,sovdt,econf,epsdgel,0,L,enert,scalpar)

         if(i.eq.inext+ibinit-nweight) then 
          do k=1,npm
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
         call convert(L,nw,nint,iconf,iconfm,ioptir)
        endif
        call reshuffhub(L,Lz,nelm2,nw,jbra,iconf
     >,table,tabpip,tabpipSz,tabpiphd,inver
     >,wsto,diagfn,diag,kelup,keldo,enert)
       endif
       if(itest.eq.2) then
        wbra=wtot
        if(iread.eq.1) then 
         call convert(L,nw,nint,iconf,iconfm,ioptir)
        endif
       endif

       if(iread.eq.1) then 
        ngg=ngg+1
        write(12) wbra,wtot,ener,ener2,diagsum,uniform_charge
     1,uniform_spin,uniform_charge_im,uniform_spin_im
     1,spin_parity(nx/2),spin_parity(nx/2-1)
        write(15) wbra,(iconfm(k),k=1,nint)
     1,(jbraw(k),k=1,nw),(wconfw(k),k=1,nw)
       else
        ngg=ngg+1
        write(12) wbra,wtot,ener,ener2,diagsum,uniform_charge
     1,uniform_spin,uniform_charge_im,uniform_spin_im
     1,spin_parity(nx/2),spin_parity(nx/2-1)
       endif
       
      enddo ! enddo for the generations

      iend=i

      if(itestr.eq.-2) then 

       write(6,*) ' Write the parameters of the final wave function '

       rewind(10)
      write(10) v,vSz,vhd,t,U,ivicr,icord,irif,nelup,neldo
     >,z,ivicn,indtn,ivicns,indtns,ivicnr,indtnr,ivicnb,indtnb
     >,vjz,vj,vjh,back,ddw,dsw,dek,dup,iesm,iesd,iesh,iback,iesfree
     >,iessw,ieskin,iesup,muf,stripe,ioptir,alphav,par,eig

      endif

      call write_seeds(nomefile)

      rewind(11)
      write(11) ngg,ngn,etry,nw,L,nbra,tbra,np,iread,iend
     >,nelup,neldo

      write(11) ((iconf(k,j),k=1,L),j=1,nw),(wconfn(j),j=1,nw)


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
c upback
c upallcorr
c upvpot
c upvpotz
c upvpothd
c uphop
c ekdk
c evaldet
c upinvhop
c upinvhop2
c upinvhopbig
c backflow
c backflowb
c countaround
c matrix
c ratiofn
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

      subroutine random(L,Lz,ztry,table,iout,indvic,gamma)
      implicit none
      integer*4 i,L,Lz,indvic,iout,indz
      real*8 try,ztry,gamma,cost
      real*8 table(*)

      try=0.d0
      i=0
      dowhile(ztry.ge.try.and.i.lt.Lz)
       i=i+1
       if(table(i).gt.0.d0) then
        try=try+table(i) 
        indz=i
       else
        cost=-gamma*table(i)
        try=try+cost
        if(cost.ne.0.d0) indz=i
       endif
      enddo

      if(i.eq.Lz) i=indz 
      indvic=(i-1)/L+1
      iout=i-(indvic-1)*L

      return
      end

      subroutine energy(Lz,table,diag,wtot,veff,ener)
      implicit none
      integer*4 Lz,i
      real*8 table(*),diag,wtot,veff,ener

      wtot=diag
      veff=0.d0
      do i=1,Lz
       wtot=wtot+table(i)
       if(table(i).lt.0.d0) veff=veff+table(i)
      enddo
      ener=-wtot

      return
      end

      subroutine initspin(nw,L,nelup,neldo,nel,iconf,kelup,keldo
     >,wconf,wpot,ipsip,psip,par,epst,iopt)
      implicit none
      integer*4 nw,L,nelup,neldo,nel,i,j,kk,ih,pos,iopt
      integer*2 iconf(L,*),kelup(L,*),keldo(L,*),keli
      integer*4 ipsip(L),info,igen,iref
      real*8 wconf(*),drand1,wpot(2*L,*),par(*)
      real*8 psip(L,L),det(2),cdet,cost,epst,dabs

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

       if(iopt.ge.0) then
        do while(ih.lt.neldo)
         pos=drand1()*L+1
c        if(iconf(pos,kk).eq.0.and.par(pos).gt.0) then
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
c        if(iconf(pos,kk).eq.0.and.par(pos).gt.0) then
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
       call dgetrf(nel,nel,psip,L,ipsip,info)
       if(info.ne.0) then
        iref=iref+1
        goto 11
       endif
       cdet=1.d0
       do i=1,nel
        cost=dabs(psip(i,i))
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

      subroutine upscratchub(L,Nw,nelm,nelup,neldo,nel
     >,tabpip,tabpipSz,tabpiphd,diag,v,vSz,vhd,U,t,iconf
     >,kelup,keldo,psip,ipsip,inver,z,iopt,ivicnb,indtnb
     >,back,iback,zet,hopupdo)
      implicit none
      integer*4 Nw,L,nelup,neldo,nel,nelm,i,j,k,h,jn,info,iopt
      integer*4 ipsip(L),indtnb(*),iback(4)
      integer*2 iconf(L,*),kelup(L,*),keldo(L,*),keli,ivicnb(L,8,*)
      real*8 U,back(L,*)
      real*8 tabpip(L,*),tabpipSz(L,*),tabpiphd(2*L,*)
     >,v(L,L),vSz(L,L),vhd(L,L),psip(nelm,2*nelm),diag(*)
     >,z(2*L,nelm),zet(nelm),t(*),hopupdo(*),inver(nelm,nelm,*)

      do k=1,nw
       if(iconf(1,k).eq.1) then
        do i=1,L
         if(iopt.ge.0) then 
          tabpip(i,k)=v(i,1)
          tabpipSz(i,k)=vSz(i,1)
          tabpiphd(i,k)=1.d0
          tabpiphd(i+L,k)=vhd(i,1)
         else
          tabpip(i,k)=v(i,1)
          tabpipSz(i,k)=vSz(i,1)
          tabpiphd(i,k)=1.d0
          tabpiphd(i+L,k)=1.d0
         endif
        enddo
       elseif(iconf(1,k).eq.-1) then
        do i=1,L
         if(iopt.ge.0) then 
          tabpip(i,k)=1.d0/v(i,1)
          tabpipSz(i,k)=vSz(i,1)
          tabpiphd(i,k)=vhd(i,1)
          tabpiphd(i+L,k)=1.d0
         else
          tabpip(i,k)=v(i,1)
          tabpipSz(i,k)=1.d0/vSz(i,1)
          tabpiphd(i,k)=1.d0
          tabpiphd(i+L,k)=1.d0
         endif
        enddo
       elseif(iconf(1,k).eq.2) then
        do i=1,L
         if(iopt.ge.0) then 
          tabpip(i,k)=1.d0
          tabpipSz(i,k)=vSz(i,1)**2
          tabpiphd(i,k)=1.d0
          tabpiphd(i+L,k)=1.d0
         else
          tabpip(i,k)=v(i,1)**2
          tabpipSz(i,k)=1.d0
          tabpiphd(i,k)=1.d0
          tabpiphd(i+L,k)=vhd(i,1)
         endif
        enddo
       elseif(iconf(1,k).eq.0) then
        do i=1,L
         if(iopt.ge.0) then 
          tabpip(i,k)=1.d0
          tabpipSz(i,k)=1.d0
          tabpiphd(i,k)=1.d0
          tabpiphd(i+L,k)=1.d0
         else
          tabpip(i,k)=1.d0
          tabpipSz(i,k)=1.d0
          tabpiphd(i,k)=vhd(i,1)
          tabpiphd(i+L,k)=1.d0
         endif
        enddo
       endif
      enddo

      do k=1,nw
       do j=2,L
        if(iconf(j,k).eq.1) then
         do i=1,L
          if(iopt.ge.0) then 
           tabpip(i,k)=tabpip(i,k)*v(i,j)
           tabpipSz(i,k)=tabpipSz(i,k)*vSz(i,j)
           tabpiphd(i+L,k)=tabpiphd(i+L,k)*vhd(i,j)
          else
           tabpip(i,k)=tabpip(i,k)*v(i,j)
           tabpipSz(i,k)=tabpipSz(i,k)*vSz(i,j)
          endif
         enddo
        elseif(iconf(j,k).eq.-1) then
         do i=1,L
          if(iopt.ge.0) then 
           tabpip(i,k)=tabpip(i,k)/v(i,j)
           tabpipSz(i,k)=tabpipSz(i,k)*vSz(i,j)
           tabpiphd(i,k)=tabpiphd(i,k)*vhd(i,j)
          else
           tabpip(i,k)=tabpip(i,k)*v(i,j)
           tabpipSz(i,k)=tabpipSz(i,k)/vSz(i,j)
          endif
         enddo
        elseif(iconf(j,k).eq.2) then 
         do i=1,L
          if(iopt.ge.0) then 
           tabpipSz(i,k)=tabpipSz(i,k)*vSz(i,j)**2
          else
           tabpip(i,k)=tabpip(i,k)*v(i,j)**2
           tabpiphd(i+L,k)=tabpiphd(i+L,k)*vhd(i,j)
          endif
         enddo
        elseif(iconf(j,k).eq.0) then 
         do i=1,L
          if(iopt.lt.0) then 
           tabpiphd(i,k)=tabpiphd(i,k)*vhd(i,j)
          endif
         enddo
        endif
       enddo
      enddo

      do k=1,nw
       do i=1,nel
        do j=1,nel
         psip(i,j)=0.d0
        enddo
       enddo

       do i=1,L

        IF (iopt.ge.0) THEN

         if(iconf(i,k).eq.2) then
          keli=kelup(i,k)
          call backflowb(L,nelm,iconf(1,k),z,zet,indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,i)
          do j=1,nel
           psip(keli,j)=zet(j)
          enddo
          keli=keldo(i,k)
          call backflowb(L,nelm,iconf(1,k),z,zet,indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,i+L)
          do j=1,nel
           psip(keli,j)=zet(j)
          enddo
         elseif(iconf(i,k).eq.1) then
          keli=kelup(i,k)
          call backflow(L,nelm,iconf(1,k),z,zet,indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,i)
          do j=1,nel
           psip(keli,j)=zet(j)
          enddo 
         elseif(iconf(i,k).eq.-1)then
          keli=keldo(i,k)
          call backflow(L,nelm,iconf(1,k),z,zet,indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,i+L)
          do j=1,nel
           psip(keli,j)=zet(j)
          enddo
         endif 

        ELSE

         if(iconf(i,k).eq.2) then
          keli=kelup(i,k)
          call backflow(L,nelm,iconf(1,k),z,zet,indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,i)
          do j=1,nel
           psip(keli,j)=zet(j)
          enddo 
          keli=keldo(i,k)
          call backflow(L,nelm,iconf(1,k),z,zet,indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,i+L)
          do j=1,nel
           psip(keli,j)=zet(j)
          enddo 
         elseif(iconf(i,k).eq.1) then
          keli=kelup(i,k)
          call backflowb(L,nelm,iconf(1,k),z,zet,indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,i)
          do j=1,nel
           psip(keli,j)=zet(j)
          enddo
         elseif(iconf(i,k).eq.-1) then
          keli=keldo(i,k)
          call backflowb(L,nelm,iconf(1,k),z,zet,indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,i+L)
          do j=1,nel
           psip(keli,j)=zet(j)
          enddo
         endif

        ENDIF

       enddo
       
       call dgetrf(nel,nel,psip,nelm,ipsip,info)
       if(info.ne.0) then
        write(6,*) ' kelup =',(kelup(i,k),i=1,L)
        write(6,*) ' keldo =',(keldo(i,k),i=1,L)
        write(6,*) 'ERROR IN DGETRF info =',info,'nw=',k
       endif

       if(info.eq.0) then 
        call dgetri(nel,psip,nelm,ipsip,psip(1,nelm+1),nelm,info)
        if(info.ne.0) write(6,*) 'ERROR IN DGETRI info =',info
       endif
       do j=1,nel
        do h=1,nel
         inver(h,j,k)=psip(h,j)
        enddo
       enddo
      enddo

      do k=1,nw
       diag(k)=0.d0
      enddo

      if(iopt.ge.0) then 
       do i=1,L
        do k=1,nw
         if(iconf(i,k).eq.1) diag(k)=diag(k)-U 
        enddo
       enddo
      else
       do i=1,L
        do k=1,nw
         if(iconf(i,k).eq.2) diag(k)=diag(k)-U 
        enddo
       enddo
      endif
 
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

      subroutine uptabtothubf(L,nel,tabpip,tabpipSz,tabpiphd
     >,i,j,iconf,vi,vj,vSzi,vSzj,vhdi,vhdj,kel,indspin,iopt)
      implicit none
      integer*4 L,i,j,k,szi,szj,nel,nelup,indspin,kels,lp,iopt
      integer*2 iconf(*),kel(L)
      real*8 tabpip(*),tabpipSz(*),tabpiphd(*)
     >,vi(*),vj(*),vSzi(*),vSzj(*),vhdi(*),vhdj(*)

      lp=L+1

      szi=iconf(i)
      szj=iconf(j)
      if(abs(iconf(i)).eq.1) szi=iconf(i)*indspin
      if(abs(iconf(j)).eq.1) szj=iconf(j)*indspin

      if(szi.eq.0.and.szj.eq.1) then
       if(iopt.ge.0) then 
        call uptabpip(L,indspin,tabpip,vi)
        call uptabpip(L,-indspin,tabpip,vj)
        call uptabpip(L,1,tabpipSz,vSzi)
        call uptabpip(L,-1,tabpipSz,vSzj)
        if(indspin.eq.1) then
         call uptabpip(L,1,tabpiphd(lp),vhdi)
         call uptabpip(L,-1,tabpiphd(lp),vhdj)
        else
         call uptabpip(L,1,tabpiphd,vhdi)
         call uptabpip(L,-1,tabpiphd,vhdj)
        endif
       else
        call uptabpip(L,1,tabpip,vi)
        call uptabpip(L,-1,tabpip,vj)
        call uptabpip(L,indspin,tabpipSz,vSzi)
        call uptabpip(L,-indspin,tabpipSz,vSzj)
        call uptabpip(L,1,tabpiphd,vhdj)
        call uptabpip(L,-1,tabpiphd,vhdi)
       endif
       kels=kel(i)
       kel(i)=kel(j)
       kel(j)=kels
       iconf(i)=iconf(j)
       iconf(j)=0
      elseif(szi.eq.0.and.szj.eq.2) then
       if(iopt.ge.0) then 
        call uptabpip(L,indspin,tabpip,vi)
        call uptabpip(L,-indspin,tabpip,vj)
        call uptabpip(L,1,tabpipSz,vSzi)
        call uptabpip(L,-1,tabpipSz,vSzj)
        if(indspin.eq.1) then
         call uptabpip(L,1,tabpiphd,vhdj)
         call uptabpip(L,1,tabpiphd(lp),vhdi)
        else
         call uptabpip(L,1,tabpiphd,vhdi)
         call uptabpip(L,1,tabpiphd(lp),vhdj)
        endif
       else
        call uptabpip(L,1,tabpip,vi)
        call uptabpip(L,-1,tabpip,vj)
        call uptabpip(L,indspin,tabpipSz,vSzi)
        call uptabpip(L,-indspin,tabpipSz,vSzj)
        call uptabpip(L,-1,tabpiphd,vhdi)
        call uptabpip(L,-1,tabpiphd(lp),vhdj)
       endif
       kels=kel(i)
       kel(i)=kel(j)
       kel(j)=kels
       iconf(i)=indspin
       iconf(j)=-iconf(i)
      elseif(szi.eq.-1.and.szj.eq.1) then
       if(iopt.ge.0) then
        call uptabpip(L,indspin,tabpip,vi)
        call uptabpip(L,-indspin,tabpip,vj)
        call uptabpip(L,1,tabpipSz,vSzi)
        call uptabpip(L,-1,tabpipSz,vSzj)
        if(indspin.eq.1) then
         call uptabpip(L,-1,tabpiphd,vhdi)
         call uptabpip(L,-1,tabpiphd(lp),vhdj)
        else
         call uptabpip(L,-1,tabpiphd,vhdj)
         call uptabpip(L,-1,tabpiphd(lp),vhdi)
        endif
       else
        call uptabpip(L,1,tabpip,vi)
        call uptabpip(L,-1,tabpip,vj)
        call uptabpip(L,indspin,tabpipSz,vSzi)
        call uptabpip(L,-indspin,tabpipSz,vSzj)
        call uptabpip(L,1,tabpiphd,vhdj)
        call uptabpip(L,1,tabpiphd(lp),vhdi)
       endif
       kels=kel(i)
       kel(i)=kel(j)
       kel(j)=kels
       iconf(i)=2
       iconf(j)=0
      elseif(szi.eq.-1.and.szj.eq.2) then
       if(iopt.ge.0) then 
        call uptabpip(L,indspin,tabpip,vi)
        call uptabpip(L,-indspin,tabpip,vj)
        call uptabpip(L,1,tabpipSz,vSzi)
        call uptabpip(L,-1,tabpipSz,vSzj)
        if(indspin.eq.1) then
         call uptabpip(L,-1,tabpiphd,vhdi)
         call uptabpip(L,1,tabpiphd,vhdj)
        else
         call uptabpip(L,-1,tabpiphd(lp),vhdi)
         call uptabpip(L,1,tabpiphd(lp),vhdj)
        endif
       else
        call uptabpip(L,1,tabpip,vi)
        call uptabpip(L,-1,tabpip,vj)
        call uptabpip(L,indspin,tabpipSz,vSzi)
        call uptabpip(L,-indspin,tabpipSz,vSzj)
        call uptabpip(L,-1,tabpiphd(lp),vhdj)
        call uptabpip(L,1,tabpiphd(lp),vhdi)
       endif
       kels=kel(i)
       kel(i)=kel(j)
       kel(j)=kels
       iconf(j)=iconf(i)
       iconf(i)=2
      endif

      return
      end

      subroutine ratiovmc(i,j,indspin,L,nelm,nel,ivic
     >,table,tabpip,tabpipSz,tabpiphd,thop,jhop
     >,z,hopupdo,zet,iconf,iconfnew,kelup,keldo
     >,iopt,izeta,iocc,sfhop,psip,ipsip,cont,cdet,matr
     >,kel,ivicnb,indtnb,iback,back)
      implicit none
      integer*4 L,i,j,jn,indspin,iopt,nelm,izeta
     >,cont,nel,k,h,info
      integer*2 ivic(L,*),ivicnb(L,8,*),iconf(L),kelup(L),keldo(L)
     >,kel(*),iconfnew(L)
      integer*4 ipsip(*),indtnb(*),iback(4)
      real*8 table,tabpip(*),tabpipSz(*),tabpiphd(*)
     >,thop(L,*),jhop(L,*),sfhop(L,*),psip(nelm,nelm)
     >,cdet(L,*),hopupdo(*),z(2*L,*),matr(L,*)
     >,zet(nelm,*),back(L,*)
      logical iocc(L)

      do k=1,L
       iconfnew(k)=iconf(k)
      enddo
     

      jn=ivic(i,j)

      table=0.d0
      cont=0

      IF(jn.eq.0) THEN
       table=0.d0
      ELSE
       if(iopt.ge.0) then 

        if(iconf(i).eq.iconf(jn)) then
         table=0.d0
        elseif(iconf(i).eq.0.and.iconf(jn).eq.1) then  ! 1
         if(indspin.eq.1) then
          iconfnew(i)=1
          iconfnew(jn)=0
          call matrix(iopt,L,nelm,nel,i,jn,jn,i,0,0
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i)/tabpiphd(jn)*thop(i,j)
         endif
         if(indspin.eq.-1) table=0.d0
        elseif(iconf(i).eq.1.and.iconf(jn).eq.0) then  ! 2
         if(indspin.eq.1) then
          iconfnew(i)=0
          iconfnew(jn)=1
          call matrix(iopt,L,nelm,nel,i,jn,i,jn,0,0
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn)/tabpiphd(i)*thop(i,j)
         endif
         if(indspin.eq.-1) table=0.d0
        elseif(iconf(i).eq.0.and.iconf(jn).eq.-1) then  ! 3
         if(indspin.eq.1) table=0.d0
         if(indspin.eq.-1) then
          iconfnew(i)=-1
          iconfnew(jn)=0
          call matrix(iopt,L,nelm,nel,i,jn,0,0,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i+L)/tabpiphd(jn+L)*thop(i,j)
         endif
        elseif(iconf(i).eq.-1.and.iconf(jn).eq.0) then  ! 4
         if(indspin.eq.1) table=0.d0
         if(indspin.eq.-1) then
          iconfnew(i)=0
          iconfnew(jn)=-1
          call matrix(iopt,L,nelm,nel,i,jn,0,0,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn+L)/tabpiphd(i+L)*thop(i,j)
         endif
        elseif(iconf(i).eq.0.and.iconf(jn).eq.2) then  ! 5
         if(indspin.eq.1) then
          iconfnew(i)=1
          iconfnew(jn)=-1
          call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn+L)*tabpiphd(i)*jhop(i,j)
         elseif(indspin.eq.-1) then
          iconfnew(i)=-1
          iconfnew(jn)=1
          call matrix(iopt,L,nelm,nel,i,jn,jn,jn,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i+L)*tabpiphd(jn)*jhop(i,j)
         elseif(indspin.eq.2) then
          iconfnew(i)=2
          iconfnew(jn)=0
          call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*(tabpipSz(i)/tabpipSz(jn))**2*sfhop(i,j)
         endif
        elseif(iconf(i).eq.2.and.iconf(jn).eq.0) then  ! 6
         if(indspin.eq.1) then
          iconfnew(i)=-1
          iconfnew(jn)=1
          call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i+L)*tabpiphd(jn)*jhop(i,j)
         elseif(indspin.eq.-1) then
          iconfnew(i)=1
          iconfnew(jn)=-1
          call matrix(iopt,L,nelm,nel,i,jn,i,i,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn+L)*tabpiphd(i)*jhop(i,j)
         elseif(indspin.eq.2) then
          iconfnew(i)=0
          iconfnew(jn)=2
          call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*(tabpipSz(jn)/tabpipSz(i))**2*sfhop(i,j)
         endif
        elseif(iconf(i).eq.1.and.iconf(jn).eq.-1) then  ! 7
         if(indspin.eq.1) then
          iconfnew(i)=0
          iconfnew(jn)=2
          call matrix(iopt,L,nelm,nel,i,jn,i,jn,jn,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >/tabpiphd(jn+L)/tabpiphd(i)*jhop(i,j)
         elseif(indspin.eq.-1) then
          iconfnew(i)=2
          iconfnew(jn)=0
          call matrix(iopt,L,nelm,nel,i,jn,i,i,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >/tabpiphd(jn+L)/tabpiphd(i)*jhop(i,j)
         endif
        elseif(iconf(i).eq.-1.and.iconf(jn).eq.1) then  ! 8
         if(indspin.eq.1) then
          iconfnew(i)=2
          iconfnew(jn)=0
          call matrix(iopt,L,nelm,nel,i,jn,jn,i,i,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >/tabpiphd(i+L)/tabpiphd(jn)*jhop(i,j)
         elseif(indspin.eq.-1) then
          iconfnew(i)=0
          iconfnew(jn)=2
          call matrix(iopt,L,nelm,nel,i,jn,jn,jn,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >/tabpiphd(i+L)/tabpiphd(jn)*jhop(i,j)
         endif
        elseif(iconf(i).eq.1.and.iconf(jn).eq.2) then  ! 9
         if(indspin.eq.1) table=0.d0
         if(indspin.eq.-1) then
          iconfnew(i)=2
          iconfnew(jn)=1
          call matrix(iopt,L,nelm,nel,i,jn,jn,jn,jn,i
     >,i,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn)/tabpiphd(i)*thop(i,j)
         endif
        elseif(iconf(i).eq.2.and.iconf(jn).eq.1) then  ! 10
         if(indspin.eq.1) table=0.d0
         if(indspin.eq.-1) then
          iconfnew(i)=1
          iconfnew(jn)=2
          call matrix(iopt,L,nelm,nel,i,jn,i,i,i,jn
     >,jn,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i)/tabpiphd(jn)*thop(i,j)
         endif
        elseif(iconf(i).eq.-1.and.iconf(jn).eq.2) then  ! 11
         if(indspin.eq.1) then
          iconfnew(i)=2
          iconfnew(jn)=-1
          call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,jn
     >,i+L,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn+L)/tabpiphd(i+L)*thop(i,j)
         endif
         if(indspin.eq.-1) table=0.d0
        elseif(iconf(i).eq.2.and.iconf(jn).eq.-1) then  ! 12
         if(indspin.eq.1) then
          iconfnew(i)=-1
          iconfnew(jn)=2
          call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,i
     >,jn+L,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i+L)/tabpiphd(jn+L)*thop(i,j)
         endif
         if(indspin.eq.-1) table=0.d0
        endif

       else

        if(iconf(i).eq.iconf(jn)) then
         table=0.d0
        elseif(iconf(i).eq.0.and.iconf(jn).eq.1) then  ! 1
         if(indspin.eq.1) then
          iconfnew(i)=1
          iconfnew(jn)=0
          call matrix(iopt,L,nelm,nel,i,jn,jn,i,0,0
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn+L)/tabpiphd(i+L)*thop(i,j)
         endif
         if(indspin.eq.-1) table=0.d0
        elseif(iconf(i).eq.1.and.iconf(jn).eq.0) then  ! 2
         if(indspin.eq.1) then
          iconfnew(i)=0
          iconfnew(jn)=1
          call matrix(iopt,L,nelm,nel,i,jn,i,jn,0,0
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i+L)/tabpiphd(jn+L)*thop(i,j)
         endif
         if(indspin.eq.-1) table=0.d0
        elseif(iconf(i).eq.0.and.iconf(jn).eq.-1) then  ! 3
         if(indspin.eq.1) table=0.d0
         if(indspin.eq.-1) then
          iconfnew(i)=-1
          iconfnew(jn)=0
          call matrix(iopt,L,nelm,nel,i,jn,0,0,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn+L)/tabpiphd(i+L)*thop(i,j)
         endif
        elseif(iconf(i).eq.-1.and.iconf(jn).eq.0) then  ! 4
         if(indspin.eq.1) table=0.d0
         if(indspin.eq.-1) then
          iconfnew(i)=0
          iconfnew(jn)=-1
          call matrix(iopt,L,nelm,nel,i,jn,0,0,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i+L)/tabpiphd(jn+L)*thop(i,j)
         endif
        elseif(iconf(i).eq.0.and.iconf(jn).eq.2) then  ! 5
         if(indspin.eq.1) then
          iconfnew(i)=1
          iconfnew(jn)=-1
          call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >/tabpiphd(i+L)/tabpiphd(jn)*jhop(i,j)
         elseif(indspin.eq.-1) then
          iconfnew(i)=-1
          iconfnew(jn)=1
          call matrix(iopt,L,nelm,nel,i,jn,jn,jn,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >/tabpiphd(i+L)/tabpiphd(jn)*jhop(i,j)
         endif
        elseif(iconf(i).eq.2.and.iconf(jn).eq.0) then  ! 6
         if(indspin.eq.1) then
          iconfnew(i)=-1
          iconfnew(jn)=1
          call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >/tabpiphd(jn+L)/tabpiphd(i)*jhop(i,j)
         elseif(indspin.eq.-1) then
          iconfnew(i)=1
          iconfnew(jn)=-1
          call matrix(iopt,L,nelm,nel,i,jn,i,i,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >/tabpiphd(jn+L)/tabpiphd(i)*jhop(i,j)
         endif
        elseif(iconf(i).eq.1.and.iconf(jn).eq.-1) then  ! 7
         if(indspin.eq.1) then
          iconfnew(i)=0
          iconfnew(jn)=2
          call matrix(iopt,L,nelm,nel,i,jn,i,jn,jn,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i+L)*tabpiphd(jn)*jhop(i,j)
         elseif(indspin.eq.-1) then
          iconfnew(i)=2
          iconfnew(jn)=0
          call matrix(iopt,L,nelm,nel,i,jn,i,i,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn+L)*tabpiphd(i)*jhop(i,j)
         elseif(indspin.eq.2) then
          iconfnew(i)=-1
          iconfnew(jn)=1
          call matrix(iopt,L,nelm,nel,i,jn,i,jn,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*(tabpipSz(jn)/tabpipSz(i))**2*sfhop(i,j)
         endif
        elseif(iconf(i).eq.-1.and.iconf(jn).eq.1) then  ! 8
         if(indspin.eq.1) then
          iconfnew(i)=2
          iconfnew(jn)=0
          call matrix(iopt,L,nelm,nel,i,jn,jn,i,i,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn+L)*tabpiphd(i)*jhop(i,j)
         elseif(indspin.eq.-1) then
          iconfnew(i)=0
          iconfnew(jn)=2
          call matrix(iopt,L,nelm,nel,i,jn,jn,jn,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i+L)*tabpiphd(jn)*jhop(i,j)
         elseif(indspin.eq.2) then
          iconfnew(i)=1
          iconfnew(jn)=-1
          call matrix(iopt,L,nelm,nel,i,jn,jn,i,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*(tabpipSz(i)/tabpipSz(jn))**2*sfhop(i,j)
         endif
        elseif(iconf(i).eq.1.and.iconf(jn).eq.2) then  ! 9
         if(indspin.eq.1) table=0.d0
         if(indspin.eq.-1) then
          iconfnew(i)=2
          iconfnew(jn)=1
          call matrix(iopt,L,nelm,nel,i,jn,jn,jn,jn,i
     >,i,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i)/tabpiphd(jn)*thop(i,j)
         endif
        elseif(iconf(i).eq.2.and.iconf(jn).eq.1) then  ! 10
         if(indspin.eq.1) table=0.d0
         if(indspin.eq.-1) then
          iconfnew(i)=1
          iconfnew(jn)=2
          call matrix(iopt,L,nelm,nel,i,jn,i,i,i,jn
     >,jn,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn)/tabpiphd(i)*thop(i,j)
         endif
        elseif(iconf(i).eq.-1.and.iconf(jn).eq.2) then  ! 11
         if(indspin.eq.1) then
          iconfnew(i)=2
          iconfnew(jn)=-1
          call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,jn
     >,i+L,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i)/tabpiphd(jn)*thop(i,j)
         endif
         if(indspin.eq.-1) table=0.d0
        elseif(iconf(i).eq.2.and.iconf(jn).eq.-1) then  ! 12
         if(indspin.eq.1) then
          iconfnew(i)=-1
          iconfnew(jn)=2
          call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,i
     >,jn+L,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
          table=table*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn)/tabpiphd(i)*thop(i,j)
         endif
         if(indspin.eq.-1) table=0.d0
        endif

       endif       

      ENDIF
 
      return
      end 

      subroutine uptablehubf(L,idim,ivc,izeta,nelm,nel,ivic
     >,table,tabpip,tabpipSz,tabpiphd,thop,jhop,iconf,iconfnew
     >,z,hopupdo,psip,ipsip,zet,kelup,keldo,kel,cdet
     >,matr,iopt,iocc,ivicnb,indtnb,iback,back)
      implicit none
      integer*4 L,idim,ivc,izeta,i,j,jn,ivcd,iopt,kj
     >,nelm,nel,cont,k,h,info
      integer*2 ivic(L,*),ivicnb(L,8,*),iconf(L),kelup(L),keldo(L)
     >,kel(*),iconfnew(L)
      integer*4 ipsip(*),indtnb(*),iback(4)
      real*8 table(L,izeta),tabpip(*),tabpipSz(*),tabpiphd(*)
     >,thop(L,*),jhop(L,*),det,z(2*L,*),zet(nelm,*)
     >,hopupdo(*),cdet(L,*),psip(nelm,nelm)
     >,matr(L,*),back(L,*)
      logical iocc(L)

      do k=1,L
       iconfnew(k)=iconf(k)
      enddo
 

      cont=0 
      ivcd=ivc*idim

      if(iopt.ge.0) then 

       do j=1,ivcd
        do i=1,L
         jn=ivic(i,j)
         IF(jn.eq.0) THEN
          table(i,j)=0.d0
          table(i,j+ivcd)=0.d0
         ELSE
          if(iconf(i).eq.iconf(jn)) then
           table(i,j)=0.d0
           table(i,j+ivcd)=0.d0
          elseif(iconf(i).eq.0.and.iconf(jn).eq.1) then  ! 1
           iconfnew(i)=1
           iconfnew(jn)=0
           call matrix(iopt,L,nelm,nel,i,jn,jn,i,0,0
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=0
           iconfnew(jn)=1
           table(i,j)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i)/tabpiphd(jn)*thop(i,j)
           table(i,j+ivcd)=0.d0
          elseif(iconf(i).eq.1.and.iconf(jn).eq.0) then  ! 2
           iconfnew(i)=0
           iconfnew(jn)=1
           call matrix(iopt,L,nelm,nel,i,jn,i,jn,0,0
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=1
           iconfnew(jn)=0
           table(i,j)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn)/tabpiphd(i)*thop(i,j)
           table(i,j+ivcd)=0.d0
          elseif(iconf(i).eq.0.and.iconf(jn).eq.-1) then  ! 3
           table(i,j)=0.d0
           iconfnew(i)=-1
           iconfnew(jn)=0
           call matrix(iopt,L,nelm,nel,i,jn,0,0,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=0
           iconfnew(jn)=-1
           table(i,j+ivcd)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i+L)/tabpiphd(jn+L)*thop(i,j)
           table(i,j+ivcd)=-table(i,j+ivcd)
          elseif(iconf(i).eq.-1.and.iconf(jn).eq.0) then  ! 4
           table(i,j)=0.d0
           iconfnew(i)=0
           iconfnew(jn)=-1
           call matrix(iopt,L,nelm,nel,i,jn,0,0,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=-1
           iconfnew(jn)=0
           table(i,j+ivcd)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn+L)/tabpiphd(i+L)*thop(i,j)
           table(i,j+ivcd)=-table(i,j+ivcd)
          elseif(iconf(i).eq.0.and.iconf(jn).eq.2) then  ! 5
           iconfnew(i)=1
           iconfnew(jn)=-1
           call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=0
           iconfnew(jn)=2
           table(i,j)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn+L)*tabpiphd(i)*jhop(i,j)
           iconfnew(i)=-1
           iconfnew(jn)=1
           call matrix(iopt,L,nelm,nel,i,jn,jn,jn,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=0
           iconfnew(jn)=2
           table(i,j+ivcd)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i+L)*tabpiphd(jn)*jhop(i,j)
           table(i,j+ivcd)=-table(i,j+ivcd)
          elseif(iconf(i).eq.2.and.iconf(jn).eq.0) then  ! 6
           iconfnew(i)=-1
           iconfnew(jn)=1
           call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=2
           iconfnew(jn)=0
           table(i,j)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i+L)*tabpiphd(jn)*jhop(i,j)
           iconfnew(i)=1
           iconfnew(jn)=-1
           call matrix(iopt,L,nelm,nel,i,jn,i,i,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=2
           iconfnew(jn)=0
           table(i,j+ivcd)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn+L)*tabpiphd(i)*jhop(i,j)
           table(i,j+ivcd)=-table(i,j+ivcd)
          elseif(iconf(i).eq.1.and.iconf(jn).eq.-1) then  ! 7
           iconfnew(i)=0
           iconfnew(jn)=2
           call matrix(iopt,L,nelm,nel,i,jn,i,jn,jn,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=1
           iconfnew(jn)=-1
           table(i,j)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >/tabpiphd(jn+L)/tabpiphd(i)*jhop(i,j)
           iconfnew(i)=2
           iconfnew(jn)=0
           call matrix(iopt,L,nelm,nel,i,jn,i,i,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=1
           iconfnew(jn)=-1
           table(i,j+ivcd)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >/tabpiphd(jn+L)/tabpiphd(i)*jhop(i,j)
           table(i,j+ivcd)=-table(i,j+ivcd)
          elseif(iconf(i).eq.-1.and.iconf(jn).eq.1) then  ! 8
           iconfnew(i)=2
           iconfnew(jn)=0
           call matrix(iopt,L,nelm,nel,i,jn,jn,i,i,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=-1
           iconfnew(jn)=1
           table(i,j)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >/tabpiphd(i+L)/tabpiphd(jn)*jhop(i,j)
           iconfnew(i)=0
           iconfnew(jn)=2
           call matrix(iopt,L,nelm,nel,i,jn,jn,jn,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=-1
           iconfnew(jn)=1
           table(i,j+ivcd)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >/tabpiphd(i+L)/tabpiphd(jn)*jhop(i,j)
           table(i,j+ivcd)=-table(i,j+ivcd)
          elseif(iconf(i).eq.1.and.iconf(jn).eq.2) then  ! 9
           table(i,j)=0.d0
           iconfnew(i)=2
           iconfnew(jn)=1
           call matrix(iopt,L,nelm,nel,i,jn,jn,jn,jn,i
     >,i,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=1
           iconfnew(jn)=2
           table(i,j+ivcd)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn)/tabpiphd(i)*thop(i,j)
           table(i,j+ivcd)=-table(i,j+ivcd)
          elseif(iconf(i).eq.2.and.iconf(jn).eq.1) then  ! 10
           table(i,j)=0.d0
           iconfnew(i)=1
           iconfnew(jn)=2
           call matrix(iopt,L,nelm,nel,i,jn,i,i,i,jn
     >,jn,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=2
           iconfnew(jn)=1
           table(i,j+ivcd)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i)/tabpiphd(jn)*thop(i,j)
           table(i,j+ivcd)=-table(i,j+ivcd)
          elseif(iconf(i).eq.-1.and.iconf(jn).eq.2) then  ! 11
           iconfnew(i)=2
           iconfnew(jn)=-1
           call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,jn
     >,i+L,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=-1
           iconfnew(jn)=2
           table(i,j)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn+L)/tabpiphd(i+L)*thop(i,j)
           table(i,j+ivcd)=0.d0
          elseif(iconf(i).eq.2.and.iconf(jn).eq.-1) then ! 12
           iconfnew(i)=-1
           iconfnew(jn)=2
           call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,i
     >,jn+L,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=2
           iconfnew(jn)=-1
           table(i,j)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i+L)/tabpiphd(jn+L)*thop(i,j)
           table(i,j+ivcd)=0.d0
          endif
         ENDIF
        enddo
       enddo

      else

       do j=1,ivcd
        do i=1,L
         jn=ivic(i,j)
         IF(jn.eq.0) THEN
          table(i,j)=0.d0
          table(i,j+ivcd)=0.d0
         ELSE
          if(iconf(i).eq.iconf(jn)) then
           table(i,j)=0.d0
           table(i,j+ivcd)=0.d0
          elseif(iconf(i).eq.0.and.iconf(jn).eq.1) then  ! 1
           iconfnew(i)=1
           iconfnew(jn)=0
           call matrix(iopt,L,nelm,nel,i,jn,jn,i,0,0
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=0
           iconfnew(jn)=1
           table(i,j)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn+L)/tabpiphd(i+L)*thop(i,j)
           table(i,j+ivcd)=0.d0
          elseif(iconf(i).eq.1.and.iconf(jn).eq.0) then  ! 2
           iconfnew(i)=0
           iconfnew(jn)=1
           call matrix(iopt,L,nelm,nel,i,jn,i,jn,0,0
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=1
           iconfnew(jn)=0
           table(i,j)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i+L)/tabpiphd(jn+L)*thop(i,j)
           table(i,j+ivcd)=0.d0
          elseif(iconf(i).eq.0.and.iconf(jn).eq.-1) then  ! 3
           table(i,j)=0.d0
           iconfnew(i)=-1
           iconfnew(jn)=0
           call matrix(iopt,L,nelm,nel,i,jn,0,0,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=0
           iconfnew(jn)=-1
           table(i,j+ivcd)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn+L)/tabpiphd(i+L)*thop(i,j)
          elseif(iconf(i).eq.-1.and.iconf(jn).eq.0) then  ! 4
           table(i,j)=0.d0
           iconfnew(i)=0
           iconfnew(jn)=-1
           call matrix(iopt,L,nelm,nel,i,jn,0,0,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=-1
           iconfnew(jn)=0
           table(i,j+ivcd)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i+L)/tabpiphd(jn+L)*thop(i,j)
          elseif(iconf(i).eq.0.and.iconf(jn).eq.2) then  ! 5
           iconfnew(i)=1
           iconfnew(jn)=-1
           call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=0
           iconfnew(jn)=2
           table(i,j)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >/tabpiphd(i+L)/tabpiphd(jn)*jhop(i,j)
           iconfnew(i)=-1
           iconfnew(jn)=1
           call matrix(iopt,L,nelm,nel,i,jn,jn,jn,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=0
           iconfnew(jn)=2
           table(i,j+ivcd)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >/tabpiphd(i+L)/tabpiphd(jn)*jhop(i,j)
          elseif(iconf(i).eq.2.and.iconf(jn).eq.0) then  ! 6
           iconfnew(i)=-1
           iconfnew(jn)=1
           call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=2
           iconfnew(jn)=0
           table(i,j)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >/tabpiphd(jn+L)/tabpiphd(i)*jhop(i,j)
           iconfnew(i)=1
           iconfnew(jn)=-1
           call matrix(iopt,L,nelm,nel,i,jn,i,i,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=2
           iconfnew(jn)=0
           table(i,j+ivcd)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >/tabpiphd(jn+L)/tabpiphd(i)*jhop(i,j)
          elseif(iconf(i).eq.1.and.iconf(jn).eq.-1) then  ! 7
           iconfnew(i)=0
           iconfnew(jn)=2
           call matrix(iopt,L,nelm,nel,i,jn,i,jn,jn,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=1
           iconfnew(jn)=-1
           table(i,j)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i+L)*tabpiphd(jn)*jhop(i,j)
           iconfnew(i)=2
           iconfnew(jn)=0
           call matrix(iopt,L,nelm,nel,i,jn,i,i,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=1
           iconfnew(jn)=-1
           table(i,j+ivcd)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn+L)*tabpiphd(i)*jhop(i,j)
          elseif(iconf(i).eq.-1.and.iconf(jn).eq.1) then  ! 8
           iconfnew(i)=2
           iconfnew(jn)=0
           call matrix(iopt,L,nelm,nel,i,jn,jn,i,i,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=-1
           iconfnew(jn)=1
           table(i,j)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn+L)*tabpiphd(i)*jhop(i,j)
           iconfnew(i)=0
           iconfnew(jn)=2
           call matrix(iopt,L,nelm,nel,i,jn,jn,jn,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=-1
           iconfnew(jn)=1
           table(i,j+ivcd)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i+L)*tabpiphd(jn)*jhop(i,j)
          elseif(iconf(i).eq.1.and.iconf(jn).eq.2) then  ! 9
           table(i,j)=0.d0
           iconfnew(i)=2
           iconfnew(jn)=1
           call matrix(iopt,L,nelm,nel,i,jn,jn,jn,jn,i
     >,i,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=1
           iconfnew(jn)=2
           table(i,j+ivcd)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(i)/tabpiphd(jn)*thop(i,j)
          elseif(iconf(i).eq.2.and.iconf(jn).eq.1) then  ! 10
           table(i,j)=0.d0
           iconfnew(i)=1
           iconfnew(jn)=2
           call matrix(iopt,L,nelm,nel,i,jn,i,i,i,jn
     >,jn,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=2
           iconfnew(jn)=1
           table(i,j+ivcd)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(jn)/tabpiphd(i)*thop(i,j)
          elseif(iconf(i).eq.-1.and.iconf(jn).eq.2) then  ! 11
           iconfnew(i)=2
           iconfnew(jn)=-1
           call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,jn
     >,i+L,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=-1
           iconfnew(jn)=2
           table(i,j)=det
     >*tabpip(i)/tabpip(jn)*tabpipSz(i)/tabpipSz(jn)
     >*tabpiphd(i)/tabpiphd(jn)*thop(i,j)
           table(i,j+ivcd)=0.d0
          elseif(iconf(i).eq.2.and.iconf(jn).eq.-1) then ! 12
           iconfnew(i)=-1
           iconfnew(jn)=2
           call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,i
     >,jn+L,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,det)
           iconfnew(i)=2
           iconfnew(jn)=-1
           table(i,j)=det
     >*tabpip(jn)/tabpip(i)*tabpipSz(jn)/tabpipSz(i)
     >*tabpiphd(jn)/tabpiphd(i)*thop(i,j)
           table(i,j+ivcd)=0.d0
          endif
         ENDIF
        enddo
       enddo

      endif

      return
      end

c====================================================================
c====================================================================
c====================================================================

      subroutine branchingo(Nw,wconfn,weight,zeta,icdiff,ipip,jbra)
      implicit none 
      integer*4 Nw,i,j,ind,icdiff,ni
      integer*4 ipip(*),jbra(Nw)
      real*8 weight,try,tryp,dble,dstep
      real*8 zeta(Nw+1),wconfn(*)

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
      do i=1,Nw
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

      subroutine reshuffhub(L,Lz,nelm2,nw,jbra,iconf
     >,table,tabpip,tabpipSz,tabpiphd,inver
     >,wsto,diagfn,diag,kelup,keldo,enert)
      implicit none
      integer*4 L,nw,j,k,ind,Lz,nelm2
      integer*4 jbra(*)
      integer*2 iconf(L,*),kelup(L,*),keldo(L,*)
      real*8 tabpip(L,*),tabpipSz(L,*),tabpiphd(2*L,*),table(Lz,*)
     >,diagfn(*),wsto(*),diag(*),inver(nelm2,*),enert(*)

      do j=1,nw
       ind=jbra(j)
       if(j.ne.ind) then
        diag(j)=diag(ind)
        diagfn(j)=diagfn(ind)
        enert(j)=enert(ind)
        wsto(j)=wsto(ind)
        call dcopy(L,tabpip(1,ind),1,tabpip(1,j),1)
        call dcopy(L,tabpipSz(1,ind),1,tabpipSz(1,j),1)
        call dcopy(2*L,tabpiphd(1,ind),1,tabpiphd(1,j),1)
        call dcopy(Lz,table(1,ind),1,table(1,j),1)
        call dcopy(nelm2,inver(1,ind),1,inver(1,j),1)
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

      subroutine reweight0(Nw,np,npm,nwm,etotw,ipip,psip,alpha,sov
     >,sovdt,econf,epst,iweight,L,enert,scalpar)
      implicit none 
      integer*4 Nw,nwm,np,npp,i,j,k,info,irank,npm,iweight,ind,L
      integer*4 ipip(*)
      real*8 dble,cost,epst
      real*8 psip(*),alpha(*),sov(npm,npm,*),sovdt(*),econf(nw,*)
     >,etotw(*),enert(*),scalpar(*)

      npp=np+1

      IF(np.ne.0) THEN 

       cost=1.d0/nw
       do i=1,npp
        do j=i,npp
         sov(i,j,1)=econf(1,i)*econf(1,j)*cost
        enddo
       enddo
       do i=2,nw
        call dsyr('U',npp,cost,econf(i,1),nw,sov,npm) 
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
         etotw(j)=etotw(j)-cost*L*enert(i)*econf(i,j)
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

      subroutine makesovop(nh,nelup,neldo,ixc,indtn,ivicn,indtns
     >,ivicns,indtnr,ivicnr,iesfree,iessw,iesup,ieskin,matrix,ss
     >,stripe,iopti,z,eig,psip,sov)
      implicit none 
      integer*4 nh,nel,nelup,neldo,k,i,j
      integer*4 ieskin,iesfree,iessw,iesup,ndim,np,stripe
     1,indtn(*),indtns(*),indtnr(*),iopti
      integer*2 ivicn(nh,8,*),ivicns(nh,8,*),ivicnr(nh,8,*)
     1,ixc(nh,*)
      real*8 ss(nh),matrix(2*nh,2*nh),psip(2*nh,*),z(2*nh,*)
     1,sov(2*nh,2*nh,*),eig(*)

      nel=nelup+neldo
      ndim=2*nh
      np=iesfree+iessw+ieskin+iesup

      do k=1,np

       if(k.le.iesfree) then 
        call upmatrix(nh,nelup,neldo,ixc,indtn,ivicn,indtns
     1,ivicns,indtnr,ivicnr,stripe,k,0,0,0,matrix,ss,iopti)

       elseif(k.le.iesfree+iessw) then 
        call upmatrix(nh,nelup,neldo,ixc,indtn,ivicn,indtns
     1,ivicns,indtnr,ivicnr,stripe,0,k-iesfree,0,0,matrix,ss,iopti)

       elseif(k.le.iesfree+iessw+iesup) then 
        call upmatrix(nh,nelup,neldo,ixc,indtn,ivicn,indtns
     1,ivicns,indtnr,ivicnr,stripe,0,0,k-iesfree-iessw,0,matrix
     1,ss,iopti)

       elseif(k.ge.iesfree+iessw+iesup+1) then 
        call upmatrix(nh,nelup,neldo,ixc,indtn,ivicn,indtns
     1,ivicns,indtnr,ivicnr,stripe,0,0,0,k-iesfree-iessw-iesup,matrix
     1,ss,iopti)

       endif

       call dgemm('N','N',ndim,ndim,ndim,1.d0,matrix,ndim,z,ndim,0.d0
     1,psip,ndim)
       call dgemm('T','N',ndim,ndim,ndim,1.d0,z,ndim,psip,ndim,0.d0
     1,matrix,ndim)

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
         matrix(j,i)=matrix(i,j)
        enddo
       enddo

       call dgemm('N','T',ndim,ndim,ndim,1.d0,matrix,ndim,z,ndim,0.d0
     1,psip,ndim)
       call dgemm('N','N',ndim,ndim,ndim,1.d0,z,ndim,psip,ndim,0.d0
     1,sov(1,1,k),ndim)

      enddo

      return
      end

      subroutine upmatrix(nh,nelup,neldo,ixc,indtn,ivicn,indtns
     >,ivicns,indtnr,ivicnr,stripe,iesfree,iessw,iesup,ieskin
     >,matrix,ss,iopti)
      implicit none 
      integer*4 nh,iq,ip,nel,nelup,neldo,k,jn
      integer*4 iesfree,iessw,iesup,ieskin,iopti,stripe
      integer*4 indtn(*),indtns(*),indtnr(*)
      integer*2 ivicn(nh,8,*),ivicns(nh,8,*),ivicnr(nh,8,*),ixc(nh,*)
      real*8 ss(nh),matrix(2*nh,2*nh),pi,dstripe

      dstripe=dfloat(stripe)

      nel=nelup+neldo

      pi=2.d0*dasin(1.d0)

      do ip=1,2*nh
       do iq=1,2*nh
        matrix(ip,iq)=0.d0
       enddo
      enddo

      if(ieskin.eq.1) then 
       do ip=1,nh
        matrix(ip,ip)=-abs(dcos(pi/dstripe*(ixc(ip,1)+0.5d0)))
        if(iopti.ge.0) then
         matrix(ip+nh,ip+nh)=abs(dcos(pi/dstripe*(ixc(ip,1)+0.5d0)))
        else
         matrix(ip+nh,ip+nh)=-abs(dcos(pi/dstripe*(ixc(ip,1)+0.5d0)))
        endif
       enddo
      endif

      if(ieskin.eq.2) then
       do ip=1,nh
        matrix(ip,ip)=-1.d0
        if(iopti.ge.0) then
         matrix(ip+nh,ip+nh)=1.d0
        else
         matrix(ip+nh,ip+nh)=-1.d0
        endif
       enddo
      endif

      if(ieskin.gt.2) then 
       do k=1,indtn(ieskin-1)
        do ip=1,nh
         jn=ivicn(ip,k,ieskin-1)
         if(jn.ne.0) then 
          if(jn.gt.0) then
           matrix(ip,jn)=matrix(ip,jn)+1.d0
           matrix(jn,ip)=matrix(jn,ip)+1.d0
           if(iopti.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)+1.d0
            matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)+1.d0
           else
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)-1.d0
            matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)-1.d0
           endif
          elseif(jn.lt.0) then
           jn=-jn
           matrix(ip,jn)=matrix(ip,jn)-1.d0
           matrix(jn,ip)=matrix(jn,ip)-1.d0
           if(iopti.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)-1.d0
            matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)-1.d0
           else
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)+1.d0
            matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)+1.d0
           endif
          endif
         endif
        enddo
       enddo
      endif

 
      if(iesfree.ne.0) then 
       do k=1,indtns(iesfree)
        do ip=1,nh
         jn=ivicns(ip,k,iesfree)
         if(jn.gt.0) then 
          matrix(ip+nh,jn)=matrix(ip+nh,jn)+1.d0
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(jn,ip+nh)=matrix(jn,ip+nh)+1.d0
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(ip,jn+nh)=matrix(ip,jn+nh)+1.d0
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(jn+nh,ip)=matrix(jn+nh,ip)+1.d0
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
         elseif(jn.lt.0) then
          jn=-jn
          matrix(ip+nh,jn)=matrix(ip+nh,jn)-1.d0
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(jn,ip+nh)=matrix(jn,ip+nh)-1.d0
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(ip,jn+nh)=matrix(ip,jn+nh)-1.d0
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(jn+nh,ip)=matrix(jn+nh,ip)-1.d0
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
         endif
        enddo
       enddo
      endif

      if(iessw.eq.1) then
       do ip=1,nh
        matrix(ip+nh,ip)=matrix(ip+nh,ip)+1.d0
        matrix(ip,ip+nh)=matrix(ip,ip+nh)+1.d0
       enddo
      endif
      if(iessw.gt.1) then 
       do k=1,indtn(iessw-1)
        do ip=1,nh
         jn=ivicn(ip,k,iessw-1)
         if(jn.gt.0) then 
          matrix(ip+nh,jn)=matrix(ip+nh,jn)+1.d0
          matrix(jn,ip+nh)=matrix(jn,ip+nh)+1.d0
          matrix(ip,jn+nh)=matrix(ip,jn+nh)+1.d0
          matrix(jn+nh,ip)=matrix(jn+nh,ip)+1.d0
         elseif(jn.lt.0) then 
          jn=-jn
          matrix(ip+nh,jn)=matrix(ip+nh,jn)-1.d0
          matrix(jn,ip+nh)=matrix(jn,ip+nh)-1.d0
          matrix(ip,jn+nh)=matrix(ip,jn+nh)-1.d0
          matrix(jn+nh,ip)=matrix(jn+nh,ip)-1.d0
         endif
        enddo
       enddo
      endif

      if(iesup.ne.0) then 
       if(iopti.ge.0) then
        do ip=1,nh
         matrix(ip,ip)=matrix(ip,ip)-ss(ip)
     <*((-1)**(ixc(ip,1)/stripe))*abs(dsin(pi/dstripe
     <*(ixc(ip,1)+0.5d0)))
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)-ss(ip)
     <*((-1)**(ixc(ip,1)/stripe))*abs(dsin(pi/dstripe
     <*(ixc(ip,1)+0.5d0)))
        enddo
       else
        do ip=1,nh
         matrix(ip+nh,ip)=matrix(ip+nh,ip)-ss(ip)
     <*((-1)**(ixc(ip,1)/stripe))*abs(dsin(pi/dstripe
     <*(ixc(ip,1)+0.5d0)))
         matrix(ip,ip+nh)=matrix(ip,ip+nh)-ss(ip)
     <*((-1)**(ixc(ip,1)/stripe))*abs(dsin(pi/dstripe
     <*(ixc(ip,1)+0.5d0)))
        enddo
       endif
      endif

      return
      end

      subroutine upback(iback,nw,L,nh2,nelm,nel,iconf,indtnb,ivicnb
     >,kelup,keldo,inver,z,hopupdo,psip,econf,iopt,iocc)
      implicit none 
      integer*4 nw,L,nh2,nelm,nel,iopt,i,j,k,jn,keli,kk,icount
      integer*4 indtnb(*),iback(4)
      integer*2 kelup(*),keldo(*),ivicnb(L,8,*),iconf(L)
      real*8 ddot
      real*8 inver(nelm,*),z(nh2,*),hopupdo(*),econf(nw,*)
     >,psip(nelm,*)
      logical iocc(L)

      icount=0

      do i=1,L
       iocc(i)=.false.
      enddo

      do kk=1,iback(1)
       icount=icount+1
       do i=1,nelm
        do j=1,nelm
         psip(i,j)=0.d0
        enddo
       enddo

       IF(iopt.ge.0) THEN

        do i=1,L
         do j=1,indtnb(kk)
          jn=ivicnb(i,j,kk)
          if((iconf(i)*iconf(jn)).eq.-1) then
           iocc(i)=.true.
           if(iconf(i).eq.1) then
            keli=kelup(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(1)*z(jn,k)
            enddo
           elseif(iconf(i).eq.-1) then
            keli=keldo(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(2)*z(jn+L,k)
            enddo
           endif
          endif
         enddo
        enddo

       ELSE

        do i=1,L
         do j=1,indtnb(kk)
          jn=ivicnb(i,j,kk)
          if((iconf(i).eq.2.and.iconf(jn).eq.0)) then
           iocc(i)=.true.
           keli=kelup(i)
           do k=1,nel
            psip(keli,k)=psip(keli,k)+hopupdo(1)*z(jn,k)
           enddo
           keli=keldo(i)
           do k=1,nel
            psip(keli,k)=psip(keli,k)+hopupdo(2)*z(jn+L,k)
           enddo
          endif
         enddo
        enddo

       ENDIF

       econf(1,icount)=0.d0

       do i=1,nelm
        econf(1,icount)=econf(1,icount)
     >  +ddot(nelm,inver(i,1),nelm,psip(1,i),1)
       enddo

      enddo

      do kk=1,iback(2)
       icount=icount+1
       do i=1,nelm
        do j=1,nelm
         psip(i,j)=0.d0
        enddo
       enddo

       IF(iopt.ge.0) THEN

        do i=1,L
         do j=1,indtnb(kk)
          jn=ivicnb(i,j,kk)
          if(iconf(i).eq.2.and.iconf(jn).eq.0) then
           keli=kelup(i)
           do k=1,nel
            psip(keli,k)=psip(keli,k)+hopupdo(1)*z(jn,k)
           enddo
           keli=keldo(i)
           do k=1,nel
            psip(keli,k)=psip(keli,k)+hopupdo(2)*z(jn+L,k)
           enddo
          endif
         enddo
        enddo

       ELSE

        do i=1,L
         do j=1,indtnb(kk)
          jn=ivicnb(i,j,kk)
          if(iconf(i)*iconf(jn).eq.-1) then
           if(iconf(i).eq.1) then
            keli=kelup(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(1)*z(jn,k)
            enddo
           elseif(iconf(i).eq.-1) then
            keli=keldo(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(2)*z(jn+L,k)
            enddo
           endif
          endif
         enddo
        enddo

       ENDIF

       econf(1,icount)=0.d0

       do i=1,nelm
        econf(1,icount)=econf(1,icount)
     >  +ddot(nelm,inver(i,1),nelm,psip(1,i),1)
       enddo

      enddo

      do kk=1,iback(3)
       icount=icount+1
       do i=1,nelm
        do j=1,nelm
         psip(i,j)=0.d0
        enddo
       enddo

       IF(iopt.ge.0) THEN

        do i=1,L
         do j=1,indtnb(kk)
          jn=ivicnb(i,j,kk)
          if(abs(iconf(i)).eq.1.and.iconf(jn).eq.0) then
           if(iconf(i).eq.1) then
            keli=kelup(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(1)*z(jn,k)
            enddo
           elseif(iconf(i).eq.-1) then
            keli=keldo(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(2)*z(jn+L,k)
            enddo
           endif
          endif
          if(iconf(i).eq.2.and.abs(iconf(jn)).eq.1) then
           if(iconf(jn).eq.-1) then
            keli=kelup(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(1)*z(jn,k)
            enddo
           elseif(iconf(jn).eq.1) then
            keli=keldo(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(2)*z(jn+L,k)
            enddo
           endif
          endif
         enddo
        enddo

       ELSE

        do i=1,L
         do j=1,indtnb(kk)
          jn=ivicnb(i,j,kk)
          if(abs(iconf(i)).eq.1.and.iconf(jn).eq.0) then
           if(iconf(i).eq.1) then
            keli=kelup(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(1)*z(jn,k)
            enddo
           elseif(iconf(i).eq.-1) then
            keli=keldo(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(2)*z(jn+L,k)
            enddo
           endif
          endif
          if(iconf(i).eq.2.and.abs(iconf(jn)).eq.1) then
           if(iconf(jn).eq.-1) then
            keli=kelup(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(1)*z(jn,k)
            enddo
           elseif(iconf(jn).eq.1) then
            keli=keldo(i)
            do k=1,nel
             psip(keli,k)=psip(keli,k)+hopupdo(2)*z(jn+L,k)
            enddo
           endif
          endif
         enddo
        enddo

       ENDIF

       econf(1,icount)=0.d0

       do i=1,nelm
        econf(1,icount)=econf(1,icount)
     >  +ddot(nelm,inver(i,1),nelm,psip(1,i),1)
       enddo

      enddo

      do kk=1,iback(4)
       icount=icount+1
       do i=1,nelm
        do j=1,nelm
         psip(i,j)=0.d0
        enddo
       enddo

       IF(iopt.ge.0) THEN

        do i=1,L
         if(iocc(i).and.iconf(i).eq.1) then
          keli=kelup(i)
          do k=1,nel
           psip(keli,k)=psip(keli,k)+z(i,k)
          enddo
         elseif(iocc(i).and.iconf(i).eq.-1) then
          keli=keldo(i)
          do k=1,nel
           psip(keli,k)=psip(keli,k)+z(i+L,k)
          enddo
         endif
         iocc(i)=.false.
        enddo

       ELSE

        do i=1,L
         if(iocc(i).and.iconf(i).eq.2) then
          keli=kelup(i)
          do k=1,nel
           psip(keli,k)=psip(keli,k)+z(i,k)
          enddo
          keli=keldo(i)
          do k=1,nel
           psip(keli,k)=psip(keli,k)+z(i+L,k)
          enddo
         endif
         iocc(i)=.false.
        enddo

       ENDIF

       econf(1,icount)=0.d0

       do i=1,nelm
        econf(1,icount)=econf(1,icount)
     >  +ddot(nelm,inver(i,1),nelm,psip(1,i),1)
       enddo

      enddo

      return
      end

      subroutine upallcorr(L,nw,nh2,nel,nelm,nelup,np
     >,sov,inver,z,amat,indtnb,ivicnb,iconf,hopupdo,back,psim
     >,kelup,keldo,econf,iopt,iback)
      implicit none 
      integer*4 L,nw,nh2,nel,nelm,nelup,np,i,j,k,h,m,jj,jn
     >,iopt,indtnb(*),iback(4),kk
      integer*2 kelup(*),keldo(*),iconf(L),ivicnb(L,8,*)
      real*8 sov(nh2,nh2,*),inver(nelm,*),econf(nw,*),amat(L,L,*)
     >,z(nh2,*),psim(nh2,nh2),hopupdo(*),back(L,*),ddot
      logical iocc1,iocc2

      do i=1,L
       do j=1,L
        amat(i,j,1)=0.d0
        amat(i,j,2)=0.d0
       enddo
      enddo

      IF(iopt.ge.0) THEN

       do i=1,L
        if(iconf(i).eq.1.or.iconf(i).eq.2) amat(i,i,1)=1.d0
        if(iconf(i).eq.-1.or.iconf(i).eq.2) amat(i,i,2)=1.d0
        iocc1=.false.
        iocc2=.false.
c       cycle over back1 + condition for back4
        do kk=1,iback(1)
         do j=1,indtnb(kk)
          jj=ivicnb(i,j,kk)
          if(iconf(i).eq.1.and.iconf(jj).eq.-1) then
           amat(i,jj,1)=back(kk,1)*hopupdo(1)
           iocc1=.true.
          endif
          if(iconf(i).eq.-1.and.iconf(jj).eq.1) then
           amat(i,jj,2)=back(kk,1)*hopupdo(2)
           iocc2=.true.
          endif
         enddo
        enddo
        if(iocc1) then
         amat(i,i,1)=amat(i,i,1)+back(1,4)-1.d0
        endif
        if(iocc2) then
         amat(i,i,2)=amat(i,i,2)+back(1,4)-1.d0
        endif
c        cycle over back2
        do kk=1,iback(2)
         do j=1,indtnb(kk)
          jj=ivicnb(i,j,kk)
          if(iconf(i).eq.2.and.iconf(jj).eq.0) then
          amat(i,jj,1)=amat(i,jj,1)+back(kk,2)*hopupdo(1)
          amat(i,jj,2)=amat(i,jj,2)+back(kk,2)*hopupdo(2)
          endif
         enddo
        enddo
c       cycle over back3
        do kk=1,iback(3)
         do j=1,indtnb(kk)
          jj=ivicnb(i,j,kk)
          if(iconf(i).eq.1.and.iconf(jj).eq.0) then
           amat(i,jj,1)=amat(i,jj,1)+back(kk,3)*hopupdo(1)
          elseif(iconf(i).eq.-1.and.iconf(jj).eq.0) then
           amat(i,jj,2)=amat(i,jj,2)+back(kk,3)*hopupdo(2)
          elseif(iconf(i).eq.2.and.iconf(jj).eq.1) then
           amat(i,jj,2)=amat(i,jj,2)+back(kk,3)*hopupdo(2)
          elseif(iconf(i).eq.2.and.iconf(jj).eq.-1) then
           amat(i,jj,1)=amat(i,jj,1)+back(kk,3)*hopupdo(1)
          endif
         enddo
        enddo
       enddo
      ELSE
       do i=1,L
        if(iconf(i).eq.1.or.iconf(i).eq.2) amat(i,i,1)=1.d0
        if(iconf(i).eq.-1.or.iconf(i).eq.2) amat(i,i,2)=1.d0
        iocc1=.false. 
c       cycle over back2
        do kk=1,iback(2)
         do j=1,indtnb(kk)
          jj=ivicnb(i,j,kk)
          if(iconf(i).eq.1.and.iconf(jj).eq.-1) then
           amat(i,jj,1)=back(kk,2)*hopupdo(1)
          elseif(iconf(i).eq.-1.and.iconf(jj).eq.1) then
           amat(i,jj,2)=back(kk,2)*hopupdo(2)
          endif
         enddo
        enddo 
c    cycle over back1+condition for back4
        do kk=1,iback(1)
         do j=1,indtnb(kk)
          jj=ivicnb(i,j,kk)
          if(iconf(i).eq.2.and.iconf(jj).eq.0) then
           amat(i,jj,1)=amat(i,jj,1)+back(kk,1)*hopupdo(1)
           amat(i,jj,2)=amat(i,jj,2)+back(kk,1)*hopupdo(2)
           iocc1=.true.
          endif
         enddo
        enddo
        if(iocc1) then
         amat(i,i,1)=amat(i,i,1)+back(1,4)-1.d0
         amat(i,i,2)=amat(i,i,2)+back(1,4)-1.d0
        endif
c     cycle over back3  
        do kk=1,iback(3)
         do j=1,indtnb(kk)
          jj=ivicnb(i,j,kk)
          if(iconf(i).eq.1.and.iconf(jj).eq.0) then
           amat(i,jj,1)=amat(i,jj,1)+back(kk,3)*hopupdo(1)
          elseif(iconf(i).eq.-1.and.iconf(jj).eq.0) then
           amat(i,jj,2)=amat(i,jj,2)+back(kk,3)*hopupdo(2)
          elseif(iconf(i).eq.2.and.iconf(jj).eq.1) then
           amat(i,jj,2)=amat(i,jj,2)+back(kk,3)*hopupdo(2)
          elseif(iconf(i).eq.2.and.iconf(jj).eq.-1) then
           amat(i,jj,1)=amat(i,jj,1)+back(kk,3)*hopupdo(1)
          endif
         enddo
        enddo
       enddo
      ENDIF
       do i=1,np
        econf(1,i)=0.d0
       enddo
       do k=1,np
        call dgemm('N','N',nh2,nh2,nh2,1.d0,sov(1,1,k),nh2,z
     >,nh2,0.d0,psim,nh2)
        do i=1,L
         do m=1,L
          if(kelup(i).ne.0) then  
           econf(1,k)=econf(1,k)
     >+amat(i,m,1)*ddot(L,inver(1,kelup(i)),1,psim(m,1),nh2)
          endif
          if(keldo(i).ne.0) then  
           econf(1,k)=econf(1,k)
     >+amat(i,m,2)*ddot(L,inver(1,keldo(i)),1,psim(m+nel,1),nh2)
          endif
         enddo
        enddo
       enddo
      return
      end

      subroutine upvpot(imis,nw,L,indtn,ivic,multirif,iconf
     >,econf,iopt)
      implicit none
      integer*4 L,i,j,k,jn,imis,nw,iopt
      integer*4 indtn(*)
      integer*2 ivic(L,8,*),iconf(*),multirif(*)
      real*8 dnn,econf(nw,*)

      if(iopt.ge.0) then 
       do k=1,imis 
        dnn=0.d0 
        if(k.eq.1) then 
         do j=1,L
          if(iconf(j).ne.2.and.iconf(j).ne.0) dnn=dnn+0.5d0
         enddo
        else
         do i=1,indtn(k-1)
          do j=1,L
           jn=abs(ivic(j,i,k-1))
           if(jn.ne.0) then
            if(iconf(j).ne.2.and.iconf(jn).ne.2) 
     >      dnn=dnn+0.5d0*multirif(k-1)*iconf(j)*iconf(jn)  
           endif
          enddo
         enddo
        endif
        econf(1,k)=dnn
       enddo
      else
       do k=1,imis 
        dnn=0.d0 
        if(k.eq.1) then
         do j=1,L
          dnn=dnn+0.5d0*iconf(j)**2
         enddo
        else
         do i=1,indtn(k-1)
          do j=1,L
           jn=abs(ivic(j,i,k-1))
           if(jn.ne.0) 
     >     dnn=dnn+0.5d0*multirif(k-1)*abs(iconf(j)*iconf(jn))
          enddo
         enddo
        endif
        econf(1,k)=dnn
       enddo
      endif

      return
      end

      subroutine upvpotz(imis,nw,L,indtn,ivic,multirif,iconf
     >,econf,iopt)
      implicit none
      integer*4 L,i,j,k,jn,imis,nw,iopt
      integer*4 indtn(*)
      integer*2 ivic(L,8,*),iconf(*),multirif(*)
      real*8 dnn,econf(nw,*)

      if(iopt.ge.0) then 
       do k=1,imis 
        dnn=0.d0 
        do i=1,indtn(k)
         do j=1,L
          jn=abs(ivic(j,i,k))
          if(jn.ne.0) 
     >    dnn=dnn+0.5d0*multirif(k)*abs(iconf(j)*iconf(jn))
         enddo
        enddo
        econf(1,k)=dnn
       enddo
      else
       do k=1,imis
        dnn=0.d0
        do i=1,indtn(k)
         do j=1,L
          jn=abs(ivic(j,i,k))
          if(jn.ne.0) then
           if(iconf(j).ne.2.and.iconf(jn).ne.2)
     >     dnn=dnn+0.5d0*multirif(k)*iconf(j)*iconf(jn)
          endif
         enddo
        enddo
        econf(1,k)=dnn
       enddo
      endif

      return
      end

      subroutine upvpothd(imis,nw,L,indtn,ivic,multirif,iconf
     >,econf,iopt)
      implicit none
      integer*4 L,i,j,k,jn,imis,nw,iopt
      integer*4 indtn(*)
      integer*2 ivic(L,8,*),iconf(*),multirif(*)
      real*8 dnn,econf(nw,*)

      if(iopt.ge.0) then 
       do k=1,imis 
        dnn=0.d0 
        do i=1,indtn(k)
         do j=1,L
          jn=abs(ivic(j,i,k))
          if(jn.ne.0) then
           if(iconf(j)*iconf(jn).eq.-1) 
     >     dnn=dnn+0.5d0*multirif(k)
          endif
         enddo
        enddo
        econf(1,k)=dnn
       enddo
      else
       do k=1,imis
        dnn=0.d0
        do i=1,indtn(k)
         do j=1,L
          jn=abs(ivic(j,i,k))
          if(jn.ne.0) then
           if((iconf(j).eq.2.and.iconf(jn).eq.0).or.
     >        (iconf(j).eq.0.and.iconf(jn).eq.2)) 
     >     dnn=dnn+0.5d0*multirif(k)
          endif
         enddo
        enddo
        econf(1,k)=dnn
       enddo
      endif

      return
      end

      subroutine uphop(L,izeth,ivic,v,vsz,vhd
     >,gh,ghsz,ghhd,t,thop,thopv,jhop,jhopv,sfhopv,ioptir)
      implicit none
      integer*4 i,k,L,jn,izeth,ioptir
      integer*2 ivic(L,*)
      real*8 thop(L,*),thopv(L,*),jhop(L,*),jhopv(L,*),sfhopv(L,*)
     >,gh(*),ghsz(*),ghhd(*),v(L,*),vsz(L,*),vhd(L,*),t(*)

      do k=1,L
       do i=1,izeth
        jn=ivic(k,i)
        if(jn.gt.0) then 
         gh(i)=v(k,k)/v(k,jn)
         ghSz(i)=vSz(k,k)/vSz(k,jn)
         ghhd(i)=vhd(k,jn)
        elseif(jn.lt.0) then
         jn=-jn
         gh(i)=-v(k,k)/v(k,jn)
         ghSz(i)=vSz(k,k)/vSz(k,jn)
         ghhd(i)=vhd(k,jn)
        else
         gh(i)=0.d0
         ghSz(i)=0.d0
        endif 
       enddo

       do i=1,izeth
        thop(k,i)=t(i)*gh(i)*ghSz(i)
        thopv(k,i)=gh(i)*ghSz(i)
        jhopv(k,i)=thopv(k,i)*ghhd(i)
        jhop(k,i)=thop(k,i)*ghhd(i)
        sfhopv(k,i)=ghSz(i)**4
       enddo
      enddo

      return
      end

      subroutine ekdk(nh,nelup,neldo,ixc,indtn,ivicn,indtns,ivicns
     >,indtnr,ivicnr,iesm,iesd,iesh,iesfree,iessw,iesup
     >,ieskin,iesslater,vjz,vj,vjh,ddw,dsw,dek,dup,psip,matrix,ss
     >,vsz,v,vhd,z,muf,stripe,ivic,izeta,t,iopti,eig)
      implicit none 
      integer*4 nh,iq,ip,nel,nelup,neldo,iopti,info,iopt,izeta
     >,i,k,j,jn,n1,n2,n3,stripe
      integer*4 iesm,iesd,iesh,iesup,iesfree,iessw,ieskin
     >,iesslater
      integer*4 indtn(*),indtns(*),indtnr(*)
      integer*2 ivicn(nh,8,*),ivicns(nh,8,*),ivicnr(nh,8,*)
     >,ixc(nh,*),ivic(nh,*)
      real*8 dup(*),ss(nh),ddw(*),dsw(*),dek(*),vj(*),vjz(*),vjh(*)
     >,psip(*),t(*),matrix(2*nh,2*nh),vsz(nh,nh),v(nh,nh),vhd(nh,nh)
     >,z(2*nh,*),eig(*),muf,dek1,pi,dstripe

      dstripe=dfloat(stripe)

      nel=nelup+neldo

      if(iesd.ne.0) then
       do i=1,nh
        do j=1,nh
         v(i,j)=1.d0
        enddo
       enddo
      endif
      if(iesm.ne.0) then
       do i=1,nh
        do j=1,nh
         vsz(i,j)=1.d0
        enddo
       enddo
      endif
      if(iesh.ne.0) then
       do i=1,nh
        do j=1,nh
         vhd(i,j)=1.d0
        enddo
       enddo
      endif

      do i=1,iesd
       if(i.eq.1) then
        do j=1,nh
         v(j,j)=dexp(vj(1))
        enddo
       else
        do k=1,indtnr(i-1)
         do j=1,nh
          jn=abs(ivicnr(j,k,i-1))
          if(jn.ne.0) then
           v(jn,j)=dexp(vj(i))
           v(j,jn)=v(jn,j)
          endif
         enddo
        enddo
       endif
      enddo

      do i=1,iesm
       do k=1,indtnr(i)
        do j=1,nh
         jn=abs(ivicnr(j,k,i))
         if(jn.ne.0) then 
          vsz(jn,j)=dexp(vjz(i))
          vsz(j,jn)=vsz(jn,j)
         endif
        enddo
       enddo
      enddo

      do i=1,iesh
       do k=1,indtnr(i)
        do j=1,nh
         jn=abs(ivicnr(j,k,i))
         if(jn.ne.0) then 
          vhd(jn,j)=dexp(vjh(i))
          vhd(j,jn)=vhd(jn,j)
         endif
        enddo
       enddo
      enddo

c Determinant

      pi=2.d0*dasin(1.d0)

      do ip=1,2*nh
       do iq=1,2*nh
        matrix(ip,iq)=0.d0
       enddo
      enddo

      if(iesslater.eq.0) return

      dek1=-1.d0
      do k=1,izeta/2
       do ip=1,nh
        jn=ivic(ip,k)
        if(jn.ne.0) then
         if(jn.gt.0) then
          matrix(ip,jn)=matrix(ip,jn)+t(k)*dek1
          matrix(jn,ip)=matrix(jn,ip)+t(k)*dek1
          if(iopti.lt.0) then
           matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)+t(k)*dek1
           matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)+t(k)*dek1
          else
           matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)-t(k)*dek1
           matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)-t(k)*dek1
          endif
         elseif(jn.lt.0) then
          jn=-jn
          matrix(ip,jn)=matrix(ip,jn)-t(k)*dek1
          matrix(jn,ip)=matrix(jn,ip)-t(k)*dek1
          if(iopti.lt.0) then
           matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)-t(k)*dek1
           matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)-t(k)*dek1
          else
           matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)+t(k)*dek1
           matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)+t(k)*dek1
          endif
         endif
        endif
       enddo
      enddo
      if(iopti.ge.0) then
       do ip=1,nh
        matrix(ip,ip)=matrix(ip,ip)-muf
        matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)+muf
       enddo
      else
       do ip=1,nh
        matrix(ip,ip)=matrix(ip,ip)-muf
        matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)-muf
       enddo
      endif

      do iq=3,ieskin
       do k=1,indtn(iq-1)
        do ip=1,nh
         jn=ivicn(ip,k,iq-1)
         if(jn.ne.0) then
          if(jn.gt.0) then 
           matrix(ip,jn)=matrix(ip,jn)+dek(iq)
           matrix(jn,ip)=matrix(jn,ip)+dek(iq)
           if(iopti.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)+dek(iq)
            matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)+dek(iq)
           else
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)-dek(iq)
            matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)-dek(iq)
           endif
          elseif(jn.lt.0) then
           jn=-jn
           matrix(ip,jn)=matrix(ip,jn)-dek(iq)
           matrix(jn,ip)=matrix(jn,ip)-dek(iq)
           if(iopti.lt.0) then
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)-dek(iq)
            matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)-dek(iq)
           else
            matrix(ip+nh,jn+nh)=matrix(ip+nh,jn+nh)+dek(iq)
            matrix(jn+nh,ip+nh)=matrix(jn+nh,ip+nh)+dek(iq)
           endif
          endif
         endif
        enddo
       enddo
      enddo

      do iq=1,iesfree
       do k=1,indtns(iq)
        do ip=1,nh
         jn=ivicns(ip,k,iq)
         if(jn.gt.0) then
          matrix(ip+nh,jn)=matrix(ip+nh,jn)+ddw(iq)
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(jn,ip+nh)=matrix(jn,ip+nh)+ddw(iq)
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(ip,jn+nh)=matrix(ip,jn+nh)+ddw(iq)
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(jn+nh,ip)=matrix(jn+nh,ip)+ddw(iq)
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
         elseif(jn.lt.0) then
          jn=-jn
          matrix(ip+nh,jn)=matrix(ip+nh,jn)-ddw(iq)
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(jn,ip+nh)=matrix(jn,ip+nh)-ddw(iq)
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(ip,jn+nh)=matrix(ip,jn+nh)-ddw(iq)
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
          matrix(jn+nh,ip)=matrix(jn+nh,ip)-ddw(iq)
     >*abs(dcos(pi/(2.d0*dstripe)*(ixc(ip,1)+1.d0+ixc(jn,1))))
         endif
        enddo
       enddo
      enddo
 

      if(iessw.gt.0) then
       do ip=1,nh
        matrix(ip+nh,ip)=matrix(ip+nh,ip)+dsw(1)
        matrix(ip,ip+nh)=matrix(ip,ip+nh)+dsw(1)
       enddo
       do iq=2,iessw
        do k=1,indtn(iq-1)
         do ip=1,nh
          jn=ivicn(ip,k,iq-1)
          if(jn.gt.0) then 
           matrix(ip+nh,jn)=matrix(ip+nh,jn)+dsw(iq)
           matrix(jn,ip+nh)=matrix(jn,ip+nh)+dsw(iq)
           matrix(ip,jn+nh)=matrix(ip,jn+nh)+dsw(iq)
           matrix(jn+nh,ip)=matrix(jn+nh,ip)+dsw(iq)
          elseif(jn.lt.0) then 
           jn=-jn
           matrix(ip+nh,jn)=matrix(ip+nh,jn)-dsw(iq)
           matrix(jn,ip+nh)=matrix(jn,ip+nh)-dsw(iq)
           matrix(ip,jn+nh)=matrix(ip,jn+nh)-dsw(iq)
           matrix(jn+nh,ip)=matrix(jn+nh,ip)-dsw(iq)
          endif
         enddo
        enddo
       enddo
      endif

      if(iesup.gt.0) then 
       if(iopti.ge.0) then 
        do ip=1,nh
         matrix(ip,ip)=matrix(ip,ip)-dup(1)*ss(ip)
     <*((-1)**(ixc(ip,1)/stripe))*abs(dsin(pi/dstripe
     <*(ixc(ip,1)+0.5d0)))
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)-dup(1)*ss(ip)
     <*((-1)**(ixc(ip,1)/stripe))*abs(dsin(pi/dstripe
     <*(ixc(ip,1)+0.5d0)))
        enddo
       else
        do ip=1,nh
         matrix(ip+nh,ip)=matrix(ip+nh,ip)-dup(1)*ss(ip)
     <*((-1)**(ixc(ip,1)/stripe))*abs(dsin(pi/dstripe
     <*(ixc(ip,1)+0.5d0)))
         matrix(ip,ip+nh)=matrix(ip,ip+nh)-dup(1)*ss(ip)
     <*((-1)**(ixc(ip,1)/stripe))*abs(dsin(pi/dstripe
     <*(ixc(ip,1)+0.5d0)))
        enddo
       endif
      endif

      if(ieskin.ge.1) then 
       if(iopti.ge.0) then 
        do ip=1,nh
         matrix(ip,ip)=matrix(ip,ip)-dek(1)
     <*abs(dcos(pi/dstripe*(ixc(ip,1)+0.5d0)))
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)+dek(1)
     <*abs(dcos(pi/dstripe*(ixc(ip,1)+0.5d0)))
        enddo
       else
        do ip=1,nh
         matrix(ip,ip)=matrix(ip,ip)-dek(1)
     <*abs(dcos(pi/dstripe*(ixc(ip,1)+0.5d0)))
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)-dek(1)
     <*abs(dcos(pi/dstripe*(ixc(ip,1)+0.5d0)))
        enddo
       endif
      endif

      if(ieskin.ge.2) then
       if(iopti.ge.0) then
        do ip=1,nh
         matrix(ip,ip)=matrix(ip,ip)-dek(2)
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)+dek(2)
        enddo
       else
        do ip=1,nh
         matrix(ip,ip)=matrix(ip,ip)-dek(2)
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)-dek(2)
        enddo
       endif
      endif


      n1=4*nh*nh+1
      n2=n1+2*nh
      n3=n2+20*nh
      call dcopy(4*nh*nh,matrix,1,psip,1)
 
      call dsyev('V','L',2*nh,psip,2*nh,psip(n1),psip(n2),20*nh,info)
 
      do i=1,2*nh
       call dcopy(2*nh,psip(2*nh*(i-1)+1),1,z(1,i),1)
       eig(i)=psip(n1+i-1)
      enddo

      return
      end

      subroutine evaldet(a,nl,nd,ipiv,det)
      implicit none
      integer*4 nl,nd,isum,i
      integer*4 ipiv(nd)
      real*8 a(nl,*),det

      det=a(1,1)
      isum=0
      if(ipiv(1).ne.1) isum=1
      do i=2,nd
       det=det*a(i,i)
       if(ipiv(i).ne.i) isum=isum+1
      enddo

      if(mod(isum,2).ne.0) det=-det

      return
      end

      subroutine upinvhop(L,nelm,nel,zet,inver,psip,epst,cdet,kel)
      implicit none
      integer*4 L,nel,i,j,nelm
      integer*2 kel(1)
      real*8 epst,zet(nelm,1),inver(nelm,*),cdet(1,1)
     >,g,ddot,psip(nelm,2)

      g=cdet(1,1)

      if(dabs(g).lt.epst)then
       return
      endif

      do j=1,nel
       psip(j,1)=inver(j,kel(1))/g
      enddo

      do j=1,nel
       psip(j,2)=ddot(nel,zet(1,1),1,inver(1,j),1)
      enddo
      psip(kel(1),2)=psip(kel(1),2)-1.d0 
        
      call dger(nel,nel,-1.d0,psip,1,psip(1,2),1,inver,nelm)

      return
      end

       subroutine upinvhop2(L,nelm,nel,zet,inver,psip,ipsip
     >,epst,cdet,kel)
       implicit none
       integer*4 L,nel,i,j,nelm,info
       integer*4 ipsip(*)
       integer*2 kel(*)
       real*8 epst,zet(nelm,2),inver(nelm,*),cdet(L,*),g,ddot
     >,psip(2*nelm,4)

       g=cdet(1,1)*cdet(2,2)-cdet(1,2)*cdet(2,1)

       if(dabs(g).lt.epst)then
        return
       endif

       g=-1.d0/g

       call dcopy(nel,inver(1,kel(1)),1,psip,1)
       call dcopy(nel,inver(1,kel(2)),1,psip(1,2),1)
       do i=1,nel
        psip(i,3)=ddot(nel,zet(1,1),1,inver(1,i),1)
       enddo
       psip(kel(1),3)=psip(kel(1),3)-1.d0
       do i=1,nel
        psip(i,4)=ddot(nel,zet(1,2),1,inver(1,i),1)
       enddo
       psip(kel(2),4)=psip(kel(2),4)-1.d0 

       do i=1,nel
        psip(i+nel,3)=cdet(2,2)*psip(i,3)-cdet(1,2)*psip(i,4)
        psip(i+nel,4)=-cdet(2,1)*psip(i,3)+cdet(1,1)*psip(i,4)
       enddo

       call dgemm('N','T',nel,nel,2,g,psip,2*nelm
     >,psip(nel+1,3),2*nelm,1.d0,inver,nelm)

       return
       end

      subroutine upinvhopbig(L,nelm,nel,zet,inver,epst,cdet,kel
     >,cont,psip,work,ipsip)
      implicit none
      integer*4 L,nel,cont,i,j,info,k,h1,h2,nelm,ind1,ind2
      integer*2 kel(*)
      integer*4 ipsip(*)
      real*8 epst,zet(nelm,cont),inver(nelm,*),cdet(L,*),g,ddot
     >,psip(2*nelm,2*cont),work(*),ginv

      call dgetrf(cont,cont,cdet,L,ipsip,info)
      if(info.ge.0) then
       call evaldet(cdet,L,cont,ipsip,g)
      else
       write(6,*)'Problem in ratiovar'
       stop
      endif

      if(dabs(g).lt.epst)then
       return
      endif
      g=-1.d0/g

      do j=1,nel
       do i=1,cont
        psip(j+nel,i+cont)=0.d0
       enddo
      enddo

      do i=1,cont
       call dcopy(nel,inver(1,kel(i)),1,psip(1,i),1)
       do j=1,nel
        psip(j,cont+i)=ddot(nel,zet(1,i),1,inver(1,j),1)
       enddo
       psip(kel(i),cont+i)=psip(kel(i),cont+i)-1.d0
      enddo
      
      call dgetri(cont,cdet,L,ipsip,work,L,info)
      if(info.ne.0)then
       write(6,*)'some problem in dgetri'
       stop
      else
       ginv=-1.d0/g
       call dgemm('N','T',nel,cont,cont,ginv,psip(1,cont+1),2*nelm
     >,cdet,L,0.d0,psip(nel+1,cont+1),2*nelm)

       call dgemm('N','T',nel,nel,cont,g,psip,2*nelm
     >,psip(nel+1,cont+1),2*nelm,1.d0,inver,nelm)
      endif

      return
      end

      subroutine backflow(L,nelm,iconf,z,zet,indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,ici)
      implicit none
      integer*4 L,nelm,ici,i,j,updo,iopt,nel,near,jn,ii
      integer*4 indtnb(*),iback(4)
      integer*2 ivicnb(L,8,*),iconf(L)
      real*8 zet(*),back(L,*),z(2*L,*),hopupdo(*),hop
      logical iocc

      updo=0
      if(ici.gt.L)then
       ici=ici-L
       updo=1
      endif

      IF(iopt.ge.0) THEN
       near=updo*2-1
       do i=1,nel
        zet(i)=0.d0
       enddo
       iocc=.false.
       do ii=1,iback(1)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         hop=hopupdo(1)
         if(updo.eq.1) hop=hopupdo(2)
         if(iconf(jn).eq.near) then
          iocc=.true.
          do i=1,nel
           zet(i)=zet(i)+back(ii,1)*hop*z(jn+updo*L,i)
          enddo
         endif
        enddo
       enddo 
       if(iocc) then
        do i=1,nel
         zet(i)=zet(i)+back(1,4)*z(ici+updo*L,i)
        enddo
       else
        do i=1,nel
         zet(i)=zet(i)+z(ici+updo*L,i)
        enddo
       endif
       do ii=1,iback(3)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         hop=hopupdo(1)
         if(updo.eq.1) hop=hopupdo(2)
         if(iconf(jn).eq.0) then
          do i=1,nel
           zet(i)=zet(i)+back(ii,3)*hop*z(jn+updo*L,i)
          enddo
         endif
        enddo
       enddo
      ELSE
       do i=1,nel
        zet(i)=0.d0
       enddo
       iocc=.false.
       do ii=1,iback(1)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         hop=hopupdo(1)
         if(updo.eq.1) hop=hopupdo(2)
         if(iconf(jn).eq.0) then
          iocc=.true.
          do i=1,nel
           zet(i)=zet(i)+back(ii,1)*hop*z(jn+updo*L,i)
          enddo
         endif
        enddo
       enddo
       if(iocc) then
        do i=1,nel
         zet(i)=zet(i)+back(1,4)*z(ici+updo*L,i)
        enddo
       else
        do i=1,nel
         zet(i)=zet(i)+z(ici+updo*L,i)
        enddo
       endif
       do ii=1,iback(3)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         hop=hopupdo(1)
         if(updo.eq.1) hop=hopupdo(2)
         if((updo.eq.0.and.iconf(jn).eq.-1).or.
     >      (updo.eq.1.and.iconf(jn).eq.1)) then
          do i=1,nel
           zet(i)=zet(i)+back(ii,3)*hop*z(jn+updo*L,i)
          enddo
         endif
        enddo
       enddo
      ENDIF

      return
      end

      subroutine backflowb(L,nelm,iconf,z,zet,indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,ici)
      implicit none
      integer*4 L,nelm,ici,i,j,updo,iopt,nel,near,jn,ii
      integer*4 indtnb(*),iback(4)
      integer*2 ivicnb(L,8,*),iconf(L)
      real*8 zet(*),back(L,*),z(2*L,*),hopupdo(*),hop

      updo=0
      if(ici.gt.L)then
       ici=ici-L
       updo=1
      endif

      IF(iopt.ge.0) THEN
       do i=1,nel
        zet(i)=z(ici+updo*L,i)
       enddo
       do ii=1,iback(2)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         hop=hopupdo(1)
         if(updo.eq.1) hop=hopupdo(2)
         if(iconf(jn).eq.0) then
          do i=1,nel
           zet(i)=zet(i)+back(ii,2)*hop*z(jn+updo*L,i)
          enddo
         endif
        enddo
       enddo 
       do ii=1,iback(3)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         hop=hopupdo(1)
         if(updo.eq.1) hop=hopupdo(2)
         if((updo.eq.0.and.iconf(jn).eq.-1).or.
     >     (updo.eq.1.and.iconf(jn).eq.1)) then
          do i=1,nel
           zet(i)=zet(i)+back(ii,3)*hop*z(jn+updo*L,i)
          enddo
         endif
        enddo
       enddo
      ELSE
       near=updo*2-1
       do i=1,nel
        zet(i)=z(ici+updo*L,i)
       enddo
       do ii=1,iback(2)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         hop=hopupdo(1)
         if(updo.eq.1) hop=hopupdo(2)
         if(iconf(jn).eq.near) then
          do i=1,nel
           zet(i)=zet(i)+back(ii,2)*hop*z(jn+updo*L,i)
          enddo
         endif
        enddo 
       enddo 
       do ii=1,iback(3)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         hop=hopupdo(1)
         if(updo.eq.1) hop=hopupdo(2)
         if(iconf(jn).eq.0) then
          do i=1,nel
           zet(i)=zet(i)+back(ii,3)*hop*z(jn+updo*L,i)
          enddo
         endif
        enddo
       enddo
      ENDIF

      return
      end

      subroutine countaround(L,nelm,iconf,indtnb,ivicnb,z,zet
     >,hopupdo,iback,back,iopt,nel,ici,icj,cont,kel,kelup,keldo
     >,iconfnew,iocc)
      implicit none
      integer*4 L,nelm,nel,iopt,i,j,ii,cont,ici,icj,k,jn,ibackmin
      integer*4 indtnb(*),iback(4)
      integer*2 ivicnb(L,8,*),iconf(L),kelup(L),keldo(L),kel(*)
     >,iconfnew(L)
      real*8 z(2*L,*),zet(nelm,*),hopupdo(*),back(L,*)
      logical iocc(L)

       do ii=1,L
        iocc(ii)=.true.
       enddo

       iocc(ici)=.false.
       iocc(icj)=.false.


c For simplicity we recalculate all the orbitals of nearest-neighbor sites
c from scratch
      IF(iopt.ge.0)THEN

       do ii=1,iback(1)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.1) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
          elseif(iconf(jn).eq.-1)then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo
       do ii=1,iback(2)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.2) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo
       ibackmin=iback(1)
       if(iback(2).lt.ibackmin) ibackmin=iback(2)
       do ii=ibackmin+1,iback(3)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.2) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          elseif(iconf(jn).eq.1) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
          elseif(iconf(jn).eq.-1)then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo

       do ii=1,iback(1)
        do j=1,indtnb(ii)
         jn=ivicnb(icj,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.1) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
          elseif(iconf(jn).eq.-1)then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo
       do ii=1,iback(2)
        do j=1,indtnb(ii)
         jn=ivicnb(icj,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.2) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo
       ibackmin=iback(1)
       if(iback(2).lt.ibackmin) ibackmin=iback(2)
       do ii=ibackmin+1,iback(3)
        do j=1,indtnb(ii)
         jn=ivicnb(icj,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.2) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          elseif(iconf(jn).eq.1) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
          elseif(iconf(jn).eq.-1)then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo

      ELSE

       do ii=1,iback(1)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.2) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo
       do ii=1,iback(2)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.1) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
          elseif(iconf(jn).eq.-1)then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo
       ibackmin=iback(1)
       if(iback(2).lt.ibackmin) ibackmin=iback(2)
       do ii=ibackmin+1,iback(3)
        do j=1,indtnb(ii)
         jn=ivicnb(ici,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.2) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          elseif(iconf(jn).eq.1) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
          elseif(iconf(jn).eq.-1)then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo

       do ii=1,iback(1)
        do j=1,indtnb(ii)
         jn=ivicnb(icj,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.2) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo
       do ii=1,iback(2)
        do j=1,indtnb(ii)
         jn=ivicnb(icj,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.1) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
          elseif(iconf(jn).eq.-1)then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo
       ibackmin=iback(1)
       if(iback(2).lt.ibackmin) ibackmin=iback(2)
       do ii=ibackmin+1,iback(3)
        do j=1,indtnb(ii)
         jn=ivicnb(icj,j,ii)
         if(iocc(jn))then
          if(iconf(jn).eq.2) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          elseif(iconf(jn).eq.1) then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=kelup(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn)
          elseif(iconf(jn).eq.-1)then
           iocc(jn)=.false.
           cont=cont+1
           kel(cont)=keldo(jn)
           call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,jn+L)
          endif 
         endif 
        enddo
       enddo

      ENDIF

      return
      end

 
      subroutine matrix(iopt,L,nelm,nel,i,jn,iupold,iupnew,idoold,idonew
     >,isite,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
      implicit none
      integer*4 iopt,L,nelm,nel,iupold,iupnew,idoold,idonew,isite
     >,cont,h,k,i,jn,izeta,info,jj
      integer*2 ivic(L,*),ivicnb(L,8,*),iconf(L),kelup(L),keldo(L)
     >,kel(*),iconfnew(L)
      integer*4 ipsip(*),indtnb(*),iback(4)
      real*8 table,cdet(L,*),hopupdo(*),z(2*L,*)
     >,matr(L,*),zet(nelm,*),psip(nelm,nelm),back(L,*),ddot
      logical iocc(L)

      cont=0

      IF(iopt.ge.0) THEN


       if(iupold.ne.0) then
        cont=cont+1
        kel(cont)=kelup(iupold)
        if(iconfnew(iupnew).eq.1) then
         call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,iupnew)
        elseif(iconfnew(iupnew).eq.2) then
         call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,iupnew)
        endif
       endif

       if(idoold.ne.0) then
        cont=cont+1
        kel(cont)=keldo(idoold)
        if(iconfnew(idonew).eq.-1) then
         call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,idonew+L)
        elseif(iconfnew(idonew).eq.2) then
         call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,idonew+L)
        endif
       endif

       if(isite.ne.0) then
        cont=cont+1
        if(isite.le.L) kel(cont)=kelup(isite)
        if(isite.gt.L) kel(cont)=keldo(isite-L)
        call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,isite)
       endif

      ELSE

       if(iupold.ne.0) then
        cont=cont+1
        kel(cont)=kelup(iupold)
        if(iconfnew(iupnew).eq.2) then
         call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,iupnew)
        elseif(iconfnew(iupnew).eq.1) then
         call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,iupnew)
        endif
       endif

       if(idoold.ne.0) then
        cont=cont+1
        kel(cont)=keldo(idoold)
        if(iconfnew(idonew).eq.2) then
         call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,idonew+L)
        elseif(iconfnew(idonew).eq.-1) then
         call backflowb(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,idonew+L)
        endif
       endif

       if(isite.ne.0) then
        cont=cont+1
        if(isite.le.L) kel(cont)=kelup(isite)
        if(isite.gt.L) kel(cont)=keldo(isite-L)
        call backflow(L,nelm,iconfnew,z,zet(1,cont),indtnb,ivicnb
     >,hopupdo,iback,back,iopt,nel,isite)
       endif

      ENDIF
 
  
      call countaround(L,nelm,iconf,indtnb,ivicnb,z,zet
     >,hopupdo,iback,back,iopt,nel,i,jn,cont,kel,kelup,keldo
     >,iconfnew,iocc)


      do h=1,cont
       do k=1,cont
        cdet(h,k)=ddot(nel,zet(1,h),1,psip(1,kel(k)),1)
        matr(h,k)=cdet(h,k)
       enddo
      enddo
      call dgetrf(cont,cont,matr,L,ipsip,info)
      if(info.ge.0) then
       call evaldet(matr,L,cont,ipsip,table)
      else
       write(6,*)'Problem in ratiovar'
       stop
      endif

      return
      end

      subroutine ratiofn(i,j,indspin,L,nelm,nel,ivic
     >,z,hopupdo,zet,iconf,iconfnew,kelup,keldo,iopt,izeta,iocc
     >,psip,ipsip,cont,cdet,matr,kel,ivicnb,indtnb,iback,back)
      implicit none
      integer*4 L,i,j,jn,indspin,iopt,nelm,izeta
     >,cont,nel,k,h,info
      integer*2 ivic(L,*),ivicnb(L,8,*),iconf(L),kelup(L),keldo(L)
     >,kel(*),iconfnew(L)
      integer*4 ipsip(*),indtnb(*),iback(4)
      real*8 psip(nelm,nelm),cdet(L,*),hopupdo(*)
     >,z(2*L,*),matr(L,*),zet(nelm,*),back(L,*),table
      logical iocc(L)

      do k=1,L
       iconfnew(k)=iconf(k)
      enddo
     
      jn=ivic(i,j)

      cont=0

      IF(jn.eq.0) THEN
       table=0.d0
      ELSE

       if(iconf(i).eq.iconf(jn)) then
        table=0.d0
       elseif(iconf(i).eq.0.and.iconf(jn).eq.1) then  ! 1
        if(indspin.eq.1) then
         iconfnew(i)=1
         iconfnew(jn)=0
         call matrix(iopt,L,nelm,nel,i,jn,jn,i,0,0
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
        if(indspin.eq.-1) table=0.d0
       elseif(iconf(i).eq.1.and.iconf(jn).eq.0) then  ! 2
        if(indspin.eq.1) then
         iconfnew(i)=0
         iconfnew(jn)=1
         call matrix(iopt,L,nelm,nel,i,jn,i,jn,0,0
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
        if(indspin.eq.-1) table=0.d0
       elseif(iconf(i).eq.0.and.iconf(jn).eq.-1) then  ! 3
        if(indspin.eq.1) table=0.d0
        if(indspin.eq.-1) then
         iconfnew(i)=-1
         iconfnew(jn)=0
         call matrix(iopt,L,nelm,nel,i,jn,0,0,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
       elseif(iconf(i).eq.-1.and.iconf(jn).eq.0) then  ! 4
        if(indspin.eq.1) table=0.d0
        if(indspin.eq.-1) then
         iconfnew(i)=0
         iconfnew(jn)=-1
         call matrix(iopt,L,nelm,nel,i,jn,0,0,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
       elseif(iconf(i).eq.0.and.iconf(jn).eq.2) then  ! 5
        if(indspin.eq.1) then
         iconfnew(i)=1
         iconfnew(jn)=-1
         call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        elseif(indspin.eq.-1) then
         iconfnew(i)=-1
         iconfnew(jn)=1
         call matrix(iopt,L,nelm,nel,i,jn,jn,jn,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
       elseif(iconf(i).eq.2.and.iconf(jn).eq.0) then  ! 6
        if(indspin.eq.1) then
         iconfnew(i)=-1
         iconfnew(jn)=1
         call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        elseif(indspin.eq.-1) then
         iconfnew(i)=1
         iconfnew(jn)=-1
         call matrix(iopt,L,nelm,nel,i,jn,i,i,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
       elseif(iconf(i).eq.1.and.iconf(jn).eq.-1) then  ! 7
        if(indspin.eq.1) then
         iconfnew(i)=0
         iconfnew(jn)=2
         call matrix(iopt,L,nelm,nel,i,jn,i,jn,jn,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        elseif(indspin.eq.-1) then
         iconfnew(i)=2
         iconfnew(jn)=0
         call matrix(iopt,L,nelm,nel,i,jn,i,i,jn,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
       elseif(iconf(i).eq.-1.and.iconf(jn).eq.1) then  ! 8
        if(indspin.eq.1) then
         iconfnew(i)=2
         iconfnew(jn)=0
         call matrix(iopt,L,nelm,nel,i,jn,jn,i,i,i
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        elseif(indspin.eq.-1) then
         iconfnew(i)=0
         iconfnew(jn)=2
         call matrix(iopt,L,nelm,nel,i,jn,jn,jn,i,jn
     >,0,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
       elseif(iconf(i).eq.1.and.iconf(jn).eq.2) then  ! 9
        if(indspin.eq.1) table=0.d0
        if(indspin.eq.-1) then
         iconfnew(i)=2
         iconfnew(jn)=1
         call matrix(iopt,L,nelm,nel,i,jn,jn,jn,jn,i
     >,i,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
       elseif(iconf(i).eq.2.and.iconf(jn).eq.1) then  ! 10
        if(indspin.eq.1) table=0.d0
        if(indspin.eq.-1) then
         iconfnew(i)=1
         iconfnew(jn)=2
         call matrix(iopt,L,nelm,nel,i,jn,i,i,i,jn
     >,jn,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
       elseif(iconf(i).eq.-1.and.iconf(jn).eq.2) then  ! 11
        if(indspin.eq.1) then
         iconfnew(i)=2
         iconfnew(jn)=-1
         call matrix(iopt,L,nelm,nel,i,jn,jn,i,jn,jn
     >,i+L,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
        if(indspin.eq.-1) table=0.d0
       elseif(iconf(i).eq.2.and.iconf(jn).eq.-1) then  ! 12
        if(indspin.eq.1) then
         iconfnew(i)=-1
         iconfnew(jn)=2
         call matrix(iopt,L,nelm,nel,i,jn,i,jn,i,i
     >,jn+L,iconf,iconfnew,kelup,keldo,z,iocc,zet,izeta,ivic
     >,hopupdo,ipsip,ivicnb,indtnb,iback,back,psip,cdet
     >,matr,cont,kel,table)
        endif
        if(indspin.eq.-1) table=0.d0
       endif

      ENDIF
 
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
