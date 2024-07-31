c---------------------------------------------------------------------
c for 2D:
c          t_x=t(1)
c          t_y=t(2)
c          t'_xy=t(3)
c          t'_yx=t(4)
c for 1D:
c          t_x=t(1)
c          t'_x=t(3)
c
c The true Deltas are twice the ones used here!
c---------------------------------------------------------------------
      program prepv
      implicit none
      integer*4 nh,nx,ny,idim,i,j,k,jn,izeta,indt,nel,nelup,neldo,ivc
      integer*4 ip,nind,nivicn,nivicnn,nivicnb,i1,ix,iy,ix0,iy0
      integer*4 iqx,iqy,iq,iesrif,iopt,ioptr,nindm,ind
     >,jj,kmax,iflagdxy,iflag,icount
     >,iesm,iesd,iesh,iesfree,iessw,iesup,ieskin,ibc,kl
     >,iback(4)
      parameter(nx=45,nh=270,idim=2,ivc=2,indt=ivc*idim,izeta=2*indt
     >,nind=nh)
      integer*2 ivic(nh,izeta),ixc(nh,idim),idiff(nh,nh)
     >,itra(nh,nh),isim(nh,8),itrasla(8*nh,nh),irif(nh)
     >,ivicns(nh,8,nind),ivicn(nh,8,nind),ivicnr(nh,8,nind)
     >,ivicnb(nh,8,nind),ivicnbs(nh,8,nind)
      integer*4 kpi(8),indtn(nind),indtns(nind),indtnr(nind)
     >,indtnb(nind),indtnbs(nind),ipip(nh),ord(9+nind)
      real*8 U,muf,rnup,rndo,xt,yt,energy
      real*8 dnrm2,dexp,dble
      real*8 v(nh,nh),vsz(nh,nh),vhd(nh,nh),alphav(indt),t(indt)
     >,delta(9),psip(4*nh*nh+22*nh)
      real*8 w(2*nh,2*nh)
      real*8 matrix(2*nh,2*nh)
      real*8 eigenval(2*nh),par(nh)
      real*8 dd(nind),dsw(nind),dup(nind),dek(nind),back(nind,4)
     >,vj(nind),vjz(nind),vjh(nind),dist(nind),distr(nind)

      namelist /symmetries/ iesrif,ioptr,iflagdxy,ibc
      namelist /couplings/ U,t
      namelist /nelectrons/ nelup,neldo
      namelist /wavefunction/ iesm,iesd,iesh,iback,iesfree,iessw,iesup
     >,ieskin,alphav,muf

      open(unit=10,file='fort.10',form='unformatted',status='unknown') 
      open(unit=11,file='tablesim',form='unformatted',status='unknown') 

      iesrif=0
      ioptr=0
      iflagdxy=0
      ibc=0
      nelup=0
      neldo=0
      iesm=0
      iesd=0
      iesh=0
      iback(1)=0
      iback(2)=0
      iback(3)=0
      iback(4)=0
      iesfree=0
      iessw=0
      iesup=0
      ieskin=0
      muf=0.d0
      do i=1,indt
       t(i)=0.d0
       alphav(i)=0.d0
      enddo
      U=0.d0
      read(5,symmetries)
      read(5,couplings)
      read(5,nelectrons)
      read(5,wavefunction)

c     iesrif=0 all the symmetryes
c     iesrif=1 no rotations
c     iesrif=2 no reflections
c     iesrif=3 no symmetries
c
c     iopt=0,-1 square lattices
c     iopt=1,-2 rectangular lattices

c     iflagdxy=0 only d_{x^2-y^2} no d_{xy}
c     iflagdxy=1 d_{x^2-y^2} on opposite sublattice and d_{xy} on the same

c     ibc=0 periodic-periodic
c     ibc=1 antiperiodic-periodic
c     ibc=2 periodic-antiperiodic
c     ibc=3 antiperiodic-antiperiodic

      if(nx.eq.nh) then
       if(iesrif.eq.0.or.iesrif.eq.2) then
        write(6,*) 'wrong iesrif',iesrif
        stop 
       endif
       if(ioptr.eq.0.or.ioptr.eq.-1) then
        write(6,*) 'wrong iopt',ioptr
        stop 
       endif
      endif

      if(ioptr.lt.0) then
       iopt=-ioptr-1
       write(6,*) 'No p.h. case'
       if(iesfree.ne.0.or.iessw.ne.0) then 
        write(6,*) 'Error no fluctuation of particles possible !!!'
        write(6,*) 'Set iesfree=iessw=0'
        stop
       endif
      else
       iopt=ioptr
       write(6,*) ' Standard p.h. case '
      endif

      do k=1,nind
       do j=1,4
        back(k,j)=0.d0
       enddo
      enddo
      back(1,4)=1.d0

      if(iesm.ne.0) then 
       read(5,*) (vjz(k),k=1,iesm)   ! Jastrow parameters Spin 
      else
       read(5,*)
      endif
      if(iesd.ne.0) then 
       read(5,*) (vj(k),k=1,iesd)    ! Jastrow parameters Density
      else
       read(5,*)
      endif
      if(iesh.ne.0) then 
       read(5,*) (vjh(k),k=1,iesh)    ! Jastrow parameters Holon-Doblon
      else
       read(5,*)
      endif
      if(iback(1).ne.0) then 
       read(5,*) (back(k,1),k=1,iback(1))
      else
       read(5,*)
      endif
      if(iback(2).ne.0) then 
       read(5,*) (back(k,2),k=1,iback(2))
      else
       read(5,*)
      endif
      if(iback(3).ne.0) then 
       read(5,*) (back(k,3),k=1,iback(3))
      else
       read(5,*)
      endif
      if(iback(4).ne.0) then 
       read(5,*) (back(k,4),k=1,iback(4))
      else
       read(5,*)
      endif
      if(iesfree.ne.0) then 
       read(5,*) (dd(k),k=1,iesfree) ! BCS d-wave parameters
      else
       read(5,*)
      endif
      if(iessw.ne.0) then 
       read(5,*) (dsw(k),k=1,iessw)  ! BCS s-wave parameters
      else
       read(5,*)
      endif
      if(iesup.ne.0) then 
       read(5,*) (dup(k),k=1,iesup)  !  AF order parameter
      else
       read(5,*)
      endif
      if(ieskin.ne.0) then
       read(5,*) (dek(k),k=1,ieskin) ! k=1: mu, k=2,ieskin: long-range hopping 
      else
       read(5,*)
      endif

      write(6,*)
      write(6,*) 'wavefunction parameters'
      write(6,*)

      write(6,*) (vjz(k),k=1,iesm)
      write(6,*) (vj(k),k=1,iesd)
      write(6,*) (vjh(k),k=1,iesh)
      write(6,*) (back(k,1),k=1,iback(1))
      write(6,*) (back(k,2),k=1,iback(2))
      write(6,*) (back(k,3),k=1,iback(3))
      write(6,*) (back(k,4),k=1,iback(4))
      write(6,*) (dd(k),k=1,iesfree)
      write(6,*) (dsw(k),k=1,iessw)
      write(6,*) (dup(k),k=1,iesup)
      write(6,*) (dek(k),k=1,ieskin)

      nel=nelup+neldo

      if(iflagdxy.eq.0) then
       do i=1,8
        kpi(i)=i
       enddo
      else ! symmetry according to d_{xy}
       kpi(1)=1
       kpi(2)=3
       kpi(3)=6
       kpi(4)=8
       kpi(5)=2
       kpi(6)=4
       kpi(7)=5
       kpi(8)=7
      endif

      call ORDINAT(IVIC,idiff,isim,irif,itra,itrasla,nx,nh,ny
     >,IXC,IXC(1,2),iesrif,iopt,ibc)
      do iq=1,nh
        par(iq)=(-1)**(ixc(iq,1)+ixc(iq,2))
      enddo

      nindm=0

      if(iopt.eq.0) then 
       write(6,*) ' Independent Momenta by symmetry ',iesrif 
       do iq=2,nh
        iqx=mod(ixc(iq,2)*nx-ixc(iq,1)*ny+10*nh-1+nh/2,nh)-nh/2+1
        iqy=mod(ixc(iq,2)*ny+ixc(iq,1)*nx+10*nh-1+nh/2,nh)-nh/2+1
        if(iesrif.ne.0.and.iqx.ge.0.and.iqy.ge.0) then       
         if(iqx.ne.0.or.iqy.ne.nh/2.or.iesrif.eq.1) then 
          nindm=nindm+1 
          write(6,*) ixc(iq,2),ixc(iq,1),'  ',iqx,iqy
         endif 
        endif 
        if(iesrif.eq.0.and.iqx.ge.0.and.iqy.ge.0.and.iqy.le.iqx) then
         nindm=nindm+1 
         write(6,*) ixc(iq,2),ixc(iq,1),'  ',iqx,iqy
        endif 
       enddo
       write(6,*) ' Independent momenta =',nindm 
      elseif(iopt.eq.1) then 
       write(6,*) ' Independent Momenta by symmetry ',iesrif 
       do iq=2,nh
        iqx=mod(ixc(iq,1)*ny+10*nh-1+nh/2,nh)-nh/2+1
        iqy=mod(ixc(iq,2)*nx+10*nh-1+nh/2,nh)-nh/2+1
        if(iesrif.ne.0.and.iqx.ge.0.and.iqy.ge.0) then       
         if(iqx.ne.0.or.iqy.ne.nh/2.or.iesrif.eq.1) then 
          nindm=nindm+1 
          write(6,*) ixc(iq,1),ixc(iq,2),'  ',iqx,iqy
         endif 
        endif 
        if(iesrif.eq.0.and.iqx.ge.0.and.iqy.ge.0.and.iqy.le.iqx) then
         nindm=nindm+1 
         write(6,*) ixc(iq,1),ixc(iq,2),'  ',iqx,iqy
        endif 
       enddo
       write(6,*) ' Independent momenta =',nindm 
      endif

      write(6,*)
      write(6,*) ' output coordinates '
      do i=1,nh
         write(6,*) i,(ixc(i,j),j=1,idim)
      enddo
        
      write(6,*)
      write(6,*) ' Table symmetries  '
      do i=1,nh
         write(6,666) i,(isim(i,j),j=1,8)
      enddo
666   format(i3,3x,8(i4,2x))

      if(ivc.le.2) then 
       if(ivc.ne.1) then 
        call prepivic(nh,nx,ivic,ivc)
       endif 

       if(ny.eq.1.and.iopt.eq.1) then 
        write(6,*) ' One-dimensional case '
        do i=1,nh
         ivic(i,2)=0
         ivic(i,4)=0
        enddo
       endif

       if(ny.eq.2.and.iopt.eq.1) then 
        write(6,*) ' Two-chain case '
        do i=1,nh
         if(ivc.eq.1) then
          if(ixc(i,2).eq.1) ivic(i,2)=0
          if(ixc(i,2).eq.0) ivic(i,4)=0
         elseif(ivc.eq.2) then
          call makeivic(nh,ixc,ivic)
         endif
        enddo
       endif
 
       write(6,*)
       write(6,*) ' Table of nearest neighbors ',ivc*idim
       do i=1,nh
        write(6,*) i,' ',(ivic(i,j),j=1,2*ivc*idim)
       enddo
      else
       write(6,*) 'Not programmed for ivc>2'
       stop
      endif

c-------------------------------------------------------------

      do i=1,nh
       ipip(i)=0
      enddo
      ipip(1)=1
 
      nivicn=0

       do i=1,nh
        if(ipip(i).eq.0) then
         nivicn=nivicn+1
         indtns(nivicn)=0
         kmax=0
         iflag=0
         do kl=1,8
          if(par(i).gt.0) then
           k=kpi(kl)
          else
           k=kl
          endif

          if(ipip(abs(isim(i,k))).eq.0) then
           kmax=kl
           ipip(abs(isim(i,k)))=1
           ipip(abs(isim(abs(isim(i,k)),3)))=1  ! mark also the reflected
           indtns(nivicn)=indtns(nivicn)+1
           if(kl.le.4) then ! the first 4 symmetries are even
            do j=1,nh
             if(isim(i,k).gt.0) then
              ivicns(j,indtns(nivicn),nivicn)=itra(j,isim(i,k))
             else
              ivicns(j,indtns(nivicn),nivicn)=-itra(j,-isim(i,k))
             endif
            enddo
           else
            do j=1,nh
             if(isim(i,k).gt.0) then
              ivicns(j,indtns(nivicn),nivicn)=-itra(j,isim(i,k))
             else
              ivicns(j,indtns(nivicn),nivicn)=itra(j,-isim(i,k))
             endif
            enddo
           endif
          else
           if(isim(i,k).eq.-i.and.k.le.4) iflag=1
          endif
         enddo
         if((kmax.gt.4.or.(iesrif.ne.0.and.iesrif.ne.2)).and.iflag.eq.0)
     >   indtns(nivicn)=-indtns(nivicn) ! dw is allowed
         do k=1,8
          ipip(abs(isim(i,k)))=1
         enddo
        endif
       enddo

      write(6,*)
      write(6,*) ' Number of independent shell =',nivicn
      write(6,*)

c table of neighbors for backflow
      do i=1,nh
       ipip(i)=0
      enddo
      ipip(1)=1
 
      nivicnb=0

       do i=1,nh
        if(ipip(i).eq.0) then
         nivicnb=nivicnb+1
         indtnbs(nivicnb)=0
         do k=1,8
          if(ipip(abs(isim(i,k))).eq.0) then
           ipip(abs(isim(i,k)))=1
           indtnbs(nivicnb)=indtnbs(nivicnb)+1
           do j=1,nh
            if(isim(i,k).gt.0) then 
             ivicnbs(j,indtnbs(nivicnb),nivicnb)=itra(j,isim(i,k))
            else
             ivicnbs(j,indtnbs(nivicnb),nivicnb)=itra(j,-isim(i,k))
            endif
           enddo
          endif
         enddo
        endif
       enddo

       if(nivicnb.ne.nivicn) then
        write(6,*) 'Something wrong with distances!'
        stop
       endif
c-------------------------------------------------------------

      if(iopt.eq.0) then 
       do i=1,nivicn
        do i1=1,9
         ord(i1) = i1
        enddo
        ip=0
        do ix=-1,1,1
         do iy=-1,1,1
          ix0=ixc(abs(ivicns(1,1,i)),1)
          iy0=ixc(abs(ivicns(1,1,i)),2)
          xt=ix0+nx*ix-ny*iy
          yt=iy0+ny*ix+nx*iy
          ip=ip+1
          delta(ip)=xt**2+yt**2
         enddo
        enddo
        call dsortx(delta,1,9,ord)
        dist(i)=dsqrt(delta(1))
       enddo
      elseif(iopt.eq.1) then 
       do i=1,nivicn
        do i1=1,9
         ord(i1) = i1
        enddo
        ip=0
        do ix=-1,1,1
         do iy=-1,1,1
          ix0=ixc(abs(ivicns(1,1,i)),1)
          iy0=ixc(abs(ivicns(1,1,i)),2)
          xt=ix0+nx*ix
          yt=iy0+ny*iy
          ip=ip+1
          delta(ip)=xt**2+yt**2
         enddo
        enddo
        call dsortx(delta,1,9,ord)
        dist(i)=dsqrt(delta(1))
       enddo
      endif

      call dsortx(dist,1,nivicn,ord)

      do i=1,nh
       do j=1,8
        do k=1,nh
         ivicn(i,j,k)=0
         ivicnr(i,j,k)=0
         ivicnb(i,j,k)=0
        enddo
       enddo
       indtn(i)=0
       indtnr(i)=0
       indtnb(i)=0
      enddo


      do i =1,nivicn
       indtn(i)=indtns(ord(i))
       do j = 1,abs(indtn(i))
        do i1 = 1,nh
         ivicn(i1,j,i) = ivicns(i1,j,ord(i))
        enddo
       enddo
      enddo

      do i =1,nivicnb
       indtnb(i)=indtnbs(ord(i))
       do j = 1,indtnb(i)
        do i1 = 1,nh
         ivicnb(i1,j,i) = ivicnbs(i1,j,ord(i))
        enddo
       enddo
      enddo


      do i=1,nh
       do j=1,8
        do k=1,nh
         ivicns(i,j,k)=0
        enddo
       enddo
       indtns(i)=0
      enddo

c-------------------------------------------------------------
c Ordering of the superconducting parameters

      write(6,*) ' Various distances: i,x,y'
      do i=1,nivicn
       write(6,777) i,dist(i),ixc(ivicn(1,1,i),1),ixc(ivicn(1,1,i),2)
      enddo
777   format(i3,2x,f14.10,2x,2(i4,2x))

      write(6,*)
      write(6,*) 'dw shell calculation'
      nivicnn=0
      do i=1,nivicn
       if(indtn(i).lt.0) then
        nivicnn=nivicnn+1
        indtns(nivicnn)=abs(indtn(i))
        do j = 1,abs(indtn(i))
         do i1 = 1,nh
          ivicns(i1,j,nivicnn) = ivicn(i1,j,i)
         enddo
        enddo
        indtn(i)=abs(indtn(i))
c it works with antiperiodic: 
c the last bonds have negative sign with respect to s-wave
        if(indtn(i).ne.1.and.iesrif.ne.3.and.iesrif.ne.1) then
         do j = indtn(i)/2+1,indtn(i)
          do i1 = 1,nh
           ivicn(i1,j,i) = -ivicn(i1,j,i)
          enddo
         enddo
        endif
        write(6,777) nivicnn,dist(i),ixc(ivicn(1,1,i),1)
     &,ixc(ivicn(1,1,i),2)
       endif
      enddo    

      if(iesfree.gt.nivicnn) then
       write(6,*) 'Too many dw parameters!',iesfree,nivicnn
       stop
      endif
      if(iessw.gt.nivicn) then
       write(6,*) 'Too many sw parameters!',iessw,nivicn
       stop
      endif

      write(6,*)
      write(6,*) 'Number of dw shells',nivicnn
      write(6,*) 

      do i=1,nivicn
       do j=1,8
        do i1=1,nh
         ivicnr(i1,j,i)=ivicn(i1,j,i)
        enddo
       enddo
       indtnr(i)=indtn(i)
       distr(i)=dist(i)
      enddo

      call ekdk(nx,ny,nh,nelup,neldo,ixc,indtn,ivicn
     >,indtns,ivicns,indtnr,ivicnr,iesm,iesd,iesh,iesfree,iessw
     >,iesup,ieskin,vjz,vj,vjh,dd,dsw,muf,dup,dek,w,eigenval,psip
     >,matrix,par,vsz,v,vhd,ivic,izeta,alphav,ioptr,ibc)

      write(6,*) '   i    Eigenvalue       Norm up         Norm down'
      do ip=1,2*nh
       rnup=dnrm2(nh,w(1,ip),1)
       if(rnup.le.1d-10) rnup=0.d0
       rndo=dnrm2(nh,w(1+nh,ip),1)
       if(rndo.le.1d-10) rndo=0.d0
       write(6,999) ip,eigenval(ip),rnup,rndo
      enddo
999   format(i5,3x,3(f14.10,2x))

      energy=0.d0
      do ip=1,nel
        energy=energy+eigenval(ip)
      enddo
      write(6,*) 
      write(6,*) ' Total GS energy    =',energy
      write(6,*) ' GS energy per site =',energy/dble(nh)

      rewind(10)
      write(10) v,vsz,vhd,t,U,ivic,ixc,irif,nelup,neldo
     >,w,ivicn,indtn,ivicns,indtns,ivicnr,indtnr,ivicnb,indtnb
     >,vjz,vj,vjh,back,dd,dsw,dek,dup,iesm,iesd,iesh,iback,iesfree
     >,iessw,ieskin,iesup,muf,ioptr,alphav,par,eigenval

      rewind(11)
      write(11) itra,isim
      close(10) 
      close(11)

      stop 
      end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine prepivic(nh,nx,ivic,ivc)
      implicit none
      integer*4 nh,nx,i,ivc 
      integer*2 ivic(nh,*),ipip1,ipip2

      if(nh.eq.nx) then
       if(ivc.eq.2) then 
        do i=1,nh
         if(ivic(i,1).gt.0) then 
          ivic(i,5)=ivic(ivic(i,1),1) ! right right 
         elseif(ivic(i,1).lt.0) then 
          ivic(i,5)=-ivic(-ivic(i,1),1) ! right right 
         else
          ivic(i,5)=0
         endif
         if(ivic(i,3).gt.0) then 
          ivic(i,7)=ivic(ivic(i,3),3) ! left left 
         elseif(ivic(i,3).lt.0) then 
          ivic(i,7)=-ivic(-ivic(i,3),3) ! left left 
         else
          ivic(i,7)=0
         endif
        enddo 
        do i=1,nh
         ipip1=ivic(i,3)
         ivic(i,3)=ivic(i,5)  ! right right
         ivic(i,4)=0
         ivic(i,5)=ipip1      ! left
         ivic(i,6)=0
         ivic(i,8)=0
        enddo
       endif
      endif

      if(nh.ne.nx) then
       if(ivc.eq.2) then 
        do i=1,nh
         if(ivic(i,1).gt.0) then 
          ivic(i,5)=ivic(ivic(i,1),2) ! right up
          ivic(i,8)=ivic(ivic(i,1),4) ! right down 
         elseif(ivic(i,1).lt.0) then 
          ivic(i,5)=-ivic(-ivic(i,1),2) ! right up
          ivic(i,8)=-ivic(-ivic(i,1),4) ! right down 
         else
          ivic(i,5)=0
          ivic(i,8)=0
         endif
         if(ivic(i,3).gt.0) then 
          ivic(i,6)=ivic(ivic(i,3),2) ! left up
          ivic(i,7)=ivic(ivic(i,3),4) ! left down 
         elseif(ivic(i,3).lt.0) then 
          ivic(i,6)=-ivic(-ivic(i,3),2) ! left up
          ivic(i,7)=-ivic(-ivic(i,3),4) ! left down 
         else
          ivic(i,6)=0
          ivic(i,7)=0
         endif
        enddo
        do i=1,nh
         ipip1=ivic(i,3)
         ipip2=ivic(i,4)
         ivic(i,3)=ivic(i,5) ! right up
         ivic(i,4)=ivic(i,6) ! left up 
         ivic(i,5)=ipip1     ! left
         ivic(i,6)=ipip2     ! down
        enddo
       endif
      endif

      return
      end

      SUBROUTINE ORDINAT(IVIC,idiff,isim,irif,itra,itrasla,nx,ndim,ly
     >,IXV,IYV,iesrif,iopt,ibc)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer*2 IXV(NDIM),IYV(NDIM),IVIC(NDIM,4),idiff(ndim,ndim) 
      integer*2 isim(ndim,8),itrasla(8*ndim,ndim),itra(ndim,ndim) 
      integer*2 irif(ndim)

      IF (iopt.eq.0) then 

c *** coordinates
       ny = nint(sqrt(float(ndim-nx*nx)))
       lx=nx
       ly=ny
       do  ii=1,ndim
        ixv(ii)=0
        iyv(ii)=0
       enddo 
       NH=NDIM
       ISHI=10*NH
       ivet=1
       do iy=0,nx+ny
        do ix=-ny,nx
         i=ix*ny-iy*nx
         j=ix*nx+iy*ny
         iop=0
         do ik=1,ndim
          i1=MOD(ixV(IK)*ny-iyV(IK)*nx-I+ISHI,NH)
          j1=MOD(ixV(IK)*nx+iyV(IK)*ny-J+ISHI,NH)
          if((I1.ne.0).OR.(J1.ne.0)) IOP=IOP+1
         enddo 
         if(iop.eq.ndim) then
          ivet=ivet+1
          ixv(ivet)=ix
          iyv(ivet)=iy
         endif
        enddo 
       enddo 
 
c *** nearest neighbors
       do ik=1,ndim
        ix=ixv(ik)
        iy=iyv(ik)
        ixp=ix+1
        iyp=iy+1
        ixm=ix-1
        iym=iy-1
        IXR=-IX
        IYR=-IY
        IXG=IY
        IYG=-IX
        i1=ixP*ny-iy*nx
        j1=ixP*nx+iy*ny
        i2=ixM*ny-iy*nx
        j2=ixM*nx+iy*ny
        i3=ix*ny-iyP*nx
        j3=ix*nx+iyP*ny
        i4=ix*ny-iyM*nx
        j4=ix*nx+iyM*ny
        i5=ixr*ny-iyr*nx
        j5=ixr*nx+iyr*ny
        i6=IXR*ny-IY*nx
        j6=IXR*nx+IY*ny
        i7=IXG*ny-IYG*nx
        j7=IXG*nx+IYG*ny
        i8=IX*ny-IYR*NX
        j8=IX*NX+IYR*NY
        do ikk=1,ndim
         IX=IXV(IKK)
         IY=IYV(IKK)
         i=ix*ny-iy*nx
         j=ix*nx+iy*ny
         ii1=i1-i+ishi
         iJ1=J1-J+ishi
         call modpar(ii1,nh,ii1s)
         call modpar(ij1,nh,ij1s)
         if((iI1.eq.0).and.(iJ1.eq.0)) then 
          ind1=ikk
          i1s=ii1s
          j1s=ij1s
         endif
         ii2=i2-i+ishi
         iJ2=J2-J+ishi
         call modpar(ii2,nh,ii2s)
         call modpar(ij2,nh,ij2s)
         if((iI2.eq.0).and.(iJ2.eq.0)) then 
          ind2=ikk
          i2s=ii2s
          j2s=ij2s
         endif
         ii3=i3-i+ishi
         iJ3=J3-J+ishi
         call modpar(ii3,nh,ii3s)
         call modpar(ij3,nh,ij3s)
         if((iI3.eq.0).and.(iJ3.eq.0)) then 
          ind3=ikk
          i3s=ii3s
          j3s=ij3s
         endif
         ii4=i4-i+ishi
         iJ4=J4-J+ishi
         call modpar(ii4,nh,ii4s)
         call modpar(ij4,nh,ij4s)
         if((iI4.eq.0).and.(iJ4.eq.0)) then 
          ind4=ikk
          i4s=ii4s
          j4s=ij4s
         endif
         ii5=i5-i+ishi
         iJ5=J5-J+ishi
         call modpar(ii5,nh,ii5s)
         call modpar(ij5,nh,ij5s)
         if((iI5.eq.0).and.(iJ5.eq.0)) then
          ind5=ikk
          i5s=ii5s
          j5s=ij5s
         endif
         ii6=i6-i+ishi
         iJ6=J6-J+ishi
         call modpar(ii6,nh,ii6s)
         call modpar(ij6,nh,ij6s)
         if((iI6.eq.0).and.(iJ6.eq.0)) then
          ind6=ikk
          i6s=ii6s
          j6s=ij6s
         endif
         ii7=i7-i+ishi
         iJ7=J7-J+ishi
         call modpar(ii7,nh,ii7s)
         call modpar(ij7,nh,ij7s)
         if((iI7.eq.0).and.(iJ7.eq.0)) then 
          ind7=ikk
          i7s=ii7s
          j7s=ij7s
         endif
         ii8=i8-i+ishi
         iJ8=J8-J+ishi
         call modpar(ii8,nh,ii8s)
         call modpar(ij8,nh,ij8s)
         if((iI8.eq.0).and.(iJ8.eq.0)) then 
          ind8=ikk
          i8s=ii8s
          j8s=ij8s
         endif
        enddo 
        if(ibc.eq.1) then 
         ivic(ik,1)=ind1*i1s
         ivic(ik,3)=ind2*i2s
         ivic(ik,2)=ind3*i3s
         ivic(ik,4)=ind4*i4s
        elseif(ibc.eq.2) then
         ivic(ik,1)=ind1*j1s
         ivic(ik,3)=ind2*j2s
         ivic(ik,2)=ind3*j3s
         ivic(ik,4)=ind4*j4s
        elseif(ibc.eq.3) then
         ivic(ik,1)=ind1*i1s*j1s
         ivic(ik,3)=ind2*i2s*j2s
         ivic(ik,2)=ind3*i3s*j3s
         ivic(ik,4)=ind4*i4s*j4s
        else
         ivic(ik,1)=ind1
         ivic(ik,3)=ind2
         ivic(ik,2)=ind3
         ivic(ik,4)=ind4
        endif
        if(ibc.eq.1) then 
         ISIM(IK,5)=IND7*i7s
         ISIM(IK,2)=IND6*i6s
         ISIM(IK,4)=IND8*i8s
         ISIM(IK,3)=IND5*i5s
        elseif(ibc.eq.2) then 
         ISIM(IK,5)=IND7*j7s
         ISIM(IK,2)=IND6*j6s
         ISIM(IK,4)=IND8*j8s
         ISIM(IK,3)=IND5*j5s
        elseif(ibc.eq.3) then 
         ISIM(IK,5)=IND7*i7s*j7s
         ISIM(IK,2)=IND6*i6s*j6s
         ISIM(IK,4)=IND8*i8s*j8s
         ISIM(IK,3)=IND5*i5s*j5s
        else
C     5=ROTATION
         ISIM(IK,5)=IND7
C     2=PARITY-X
         ISIM(IK,2)=IND6
C     4=PARITY-Y
         ISIM(IK,4)=IND8
C     3=inversion 
         ISIM(IK,3)=IND5
        endif
       enddo 

       if(iesrif.eq.2.or.iesrif.eq.3) then 
        do i=1,ndim
         isim(i,2)=i
         isim(i,4)=i
        enddo
       endif
       if(iesrif.eq.1.or.iesrif.eq.3) then 
        do i=1,ndim
         isim(i,5)=i
        enddo
       endif
 
c *** table of idiff
       do i=1,ndim
        do j=1,ndim
         ixd=ixv(i)-ixv(j)
         iyd=iyv(i)-iyv(j)
         id=ixd*ny-iyd*nx
         jd=ixd*nx+iyd*ny
         do ikk=1,ndim
          IX=IXV(IKK)
          IY=IYV(IKK)
          i1=ix*ny-iy*nx
          j1=ix*nx+iy*ny
          iI1=mod(i1-Id+ishi,nh)
          iJ1=mod(J1-Jd+ishi,nh)
          if((iI1.eq.0).and.(iJ1.eq.0)) ind1=ikk
         enddo
         idiff(i,j)=ind1
        enddo
       enddo
 
c *** table of translations
       do ik=1,ndim
        ix=ixv(ik)
        iy=iyv(ik)
        do ikt=1,ndim
         ixt=ixv(ikt)+ix
         iyt=iyv(ikt)+iy
         i1=IXt*ny-IYt*nx
         j1=IXt*nx+IYt*ny
         do ikk=1,ndim
          IXxx=IXV(IKK)
          IYyy=IYV(IKK)
          i=ixxx*ny-iyyy*nx-I1+ISHI
          j=ixxx*nx+iyyy*ny-J1+ISHI
          call modpar(i,nh,iis)
          call modpar(j,nh,ijs)
          if((I.eq.0).and.(J.eq.0)) then 
           IND1=IKK
           is=iis
           js=ijs
          endif
         enddo
         if(ibc.eq.1) then 
          itra(ik,ikt)=ind1*is
         elseif(ibc.eq.2) then 
          itra(ik,ikt)=ind1*js
         elseif(ibc.eq.3) then 
          itra(ik,ikt)=ind1*is*js
         else
          itra(ik,ikt)=ind1
         endif
        enddo
       enddo
 
      ELSE
 
c *** coordinates
c *** nearest neighbors
c *** symmetric points
c *** table of translations
       ny=ndim/nx
       lx=nx
       ly=ny
       do j=0,ly-1
        do i=0,lx-1
         ind=j*lx+i+1
         ixv(ind)=i
         iyv(ind)=j
        enddo
       enddo
       do i=1,ndim
        iy=(i-1)/lx
        ix=i-iy*lx-1
        do j=1,ndim
         iytp=(j-1)/lx
         ixtp=j-iytp*lx-1
         ixt=mod(ix+ixtp,lx)
         iyt=mod(iy+iytp,ly)
         ind=iyt*lx+ixt+1
         if(ibc.eq.1) then 
          if(ixt.ne.ix+ixtp) ind=-ind
         endif
         if(ibc.eq.2) then 
          if(iyt.ne.iy+iytp) ind=-ind
         endif
         if(ibc.eq.3) then 
          if(ixt.ne.ix+ixtp) ind=-ind
          if(iyt.ne.iy+iytp) ind=-ind
         endif
         itra(i,j)=ind
        enddo
       enddo
       do i=1,ndim
        iy=(i-1)/lx
        ix=i-iy*lx-1
        ixt=mod(lx-ix,lx)
        indpx=lx*iy+ixt+1
        if((ibc.eq.1.or.ibc.eq.3).and.ix.ne.0) indpx=-indpx
        iyt=mod(ly-iy,ly)
        indpy=lx*iyt+ix+1
        if((ibc.eq.2.or.ibc.eq.3).and.iy.ne.0) indpy=-indpy
        ixr=mod(lx-iy,lx)
        iyr=ix
        indr=iyr*lx+ixr+1
        if((ibc.eq.2.or.ibc.eq.3).and.iy.ne.0) indr=-indr
        isim(i,2)=indpx
        isim(i,4)=indpy
        isim(i,5)=indr
       enddo
       do i=1,ndim
        if(isim(i,2).gt.0) then
         isim(i,3)=isim(isim(i,2),4)
        else
         isim(i,3)=-isim(-isim(i,2),4)
        endif
       enddo
       if(iesrif.eq.2.or.iesrif.eq.3) then
        do i=1,ndim
         isim(i,2)=i
         isim(i,4)=i
        enddo
       endif
       if(iesrif.eq.1.or.iesrif.eq.3) then 
        do i=1,ndim
         isim(i,5)=i
        enddo
       endif

c *** nearest neighbors
       do i=1,ndim
        itest=i/lx
        if(itest*lx.eq.i)then
         ir=lx
        else
         ir=0
        endif
        if(itest*lx.eq.i-1)then
         il=lx
        else
         il=0
        endif
        if(i.gt.ndim-lx)then
         iup=ndim
        else
         iup=0
        endif
        if(i.le.lx)then
         idow=ndim
        else
         idow=0
        endif
        ivic(i,1)=i+1-ir        ! right
        ivic(i,3)=i-1+il        ! left
        ivic(i,2)=i+lx-iup      ! up
        ivic(i,4)=i-lx+idow     ! down
        if(ibc.eq.1) then
         if(ir.ne.0) ivic(i,1)=-ivic(i,1)
         if(il.ne.0) ivic(i,3)=-ivic(i,3)
        endif
        if(ibc.eq.2) then
         if(iup.ne.0) ivic(i,2)=-ivic(i,2)
         if(idow.ne.0) ivic(i,4)=-ivic(i,4)
        endif
        if(ibc.eq.3) then
         if(ir.ne.0) ivic(i,1)=-ivic(i,1)
         if(il.ne.0) ivic(i,3)=-ivic(i,3)
         if(iup.ne.0) ivic(i,2)=-ivic(i,2)
         if(idow.ne.0) ivic(i,4)=-ivic(i,4)
        endif
       enddo    

c *** table of idiff
       do i=1,ndim
        iy=(i-1)/lx
        ix=i-iy*lx-1
        do j=1,ndim
         iyt=(j-1)/lx
         ixt=j-iyt*lx-1
         ixt=mod(ix-ixt+lx,lx)
         iyt=mod(iy-iyt+ly,ly)
         ind=iyt*lx+ixt+1
         idiff(i,j)=ind
        enddo
       enddo

      ENDIF

      do I=1,NDIM
       if(isim(i,5).gt.0) then 
        ISIM(I,6)=ISIM(ISIM(I,5),2)
       else
        ISIM(I,6)=-ISIM(-ISIM(I,5),2)
       endif
       if(iesrif.eq.1.or.iesrif.eq.3) isim(i,6)=i
      enddo
      do I=1,NDIM
       if(isim(i,3).gt.0) then 
        ISIM(I,7)=ISIM(ISIM(I,3),5)
       else
        ISIM(I,7)=-ISIM(-ISIM(I,3),5)
       endif
       if(iesrif.eq.1.or.iesrif.eq.3) isim(i,7)=i
      enddo
      do I=1,NDIM
       if(isim(i,7).gt.0) then 
        ISIM(I,8)=ISIM(ISIM(I,7),2)
       else
        ISIM(I,8)=-ISIM(-ISIM(I,7),2)
       endif
       if(iesrif.eq.1.or.iesrif.eq.3) isim(i,8)=i
      enddo
      do I=1,NDIM
       ISIM(I,1)=I
      enddo
      
      isp=8
      
      DO ISPA=1,isp
       DO ITR=1,NDIM
        ITOT=NDIM*(ISPA-1)+ITR
        DO ISITE=1,NDIM
         if(itra(isite,itr).gt.0) then 
          ITRASLA(ITOT,ISITE)=ISIM(ITRA(ISITE,ITR),ISPA)
         else
          ITRASLA(ITOT,ISITE)=-ISIM(-ITRA(ISITE,ITR),ISPA)
         endif
        ENDDO
       ENDDO
      ENDDO

      do i=1,ndim
       irif(i)=isim(i,3)
      enddo

      return
      end

      subroutine modpar(i,L,ipar)
      implicit none
      integer*4 i,L,ipar

      ipar=i/L
      i=i-ipar*L

      if(mod(ipar,2).eq.0) then
       ipar=1
      else
       ipar=-1
      endif

      return
      end

      subroutine ekdk(nx,ny,nh,nelup,neldo,ixc,indtn,ivicn
     >,indtns,ivicns,indtnr,ivicnr,iesm,iesd,iesh,iesfree,iessw
     >,iesup,ieskin,vjz,vj,vjh,dd,dsw,muf,dup,dek,w,eigenval,psip
     >,matrix,ss,vsz,v,vhd,ivic,izeta,t,iopt,ibc)
      implicit none 
      integer*4 nh,iq,ip,iqx,iqy,nx,ny,nel,nelup,neldo,izeta
     >,i,k,j,jn,ind,iopt,ibc,n1,n2,n3
      integer*4 icount,info,iesm,iesd,iesh
     >,iesfree,iessw,iesup,ieskin,indtn(*),indtns(*),indtnr(*)
      integer*2 ivicn(nh,8,*),ivicns(nh,8,*),ivicnr(nh,8,*)
     >,ixc(nh,*),ivic(nh,*)
      real*8 dek(*),dup(*),dd(*),dsw(*),w(2*nh,2*nh)
     >,ss(nh),psip(*),vj(*),vjz(*),vjh(*),dek1
     >,muf,matrix(2*nh,2*nh),eigenval(2*nh),t(*)
     >,vsz(nh,nh),v(nh,nh),vhd(nh,nh)

c iopt>=0  with p.h.
c iopt<0  with no p.h.

c Prepare the Jastrow terms

       do i=1,nh
        do j=1,nh
         vsz(i,j)=1.d0
         v(i,j)=1.d0
         vhd(i,j)=1.d0
        enddo
       enddo
       
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

c Prepare the Slater determinant
      
      write(6,*)
      write(6,*) ' Chosen mu       =',muf
      if(ieskin.ge.1) then
       write(6,*) ' Renormalized mu =',muf+dek(1)
      endif
      write(6,*)
         
      do ip=1,2*nh
       do iq=1,2*nh
        matrix(ip,iq)=0.d0
       enddo
      enddo

      dek1=-1.d0
      do k=1,izeta/2
       do ip=1,nh
        jn=ivic(ip,k)
        if(jn.ne.0) then
         if(jn.gt.0) then 
          matrix(ip,jn)=matrix(ip,jn)+t(k)*dek1
          matrix(jn,ip)=matrix(jn,ip)+t(k)*dek1
          if(iopt.lt.0) then
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
          if(iopt.lt.0) then
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

      do iq=2,ieskin
       do k=1,indtn(iq)
        do ip=1,nh
         jn=ivicn(ip,k,iq)
         if(jn.ne.0) then
          if(jn.gt.0) then
           matrix(ip,jn)=matrix(ip,jn)+dek(iq)
           matrix(jn,ip)=matrix(jn,ip)+dek(iq)
           if(iopt.lt.0) then
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
           if(iopt.lt.0) then
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
          matrix(ip+nh,jn)=matrix(ip+nh,jn)+dd(iq)
          matrix(jn,ip+nh)=matrix(jn,ip+nh)+dd(iq)
          matrix(ip,jn+nh)=matrix(ip,jn+nh)+dd(iq)
          matrix(jn+nh,ip)=matrix(jn+nh,ip)+dd(iq)
         elseif(jn.lt.0) then 
          jn=-jn
          matrix(ip+nh,jn)=matrix(ip+nh,jn)-dd(iq)
          matrix(jn,ip+nh)=matrix(jn,ip+nh)-dd(iq)
          matrix(ip,jn+nh)=matrix(ip,jn+nh)-dd(iq)
          matrix(jn+nh,ip)=matrix(jn+nh,ip)-dd(iq)
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

      if(iopt.ge.0) then 
       if(iesup.gt.0) then  ! magnetization along z
        do ip=1,nh
         matrix(ip,ip)=matrix(ip,ip)-dup(1)*ss(ip)
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)-dup(1)*ss(ip)
        enddo
       endif
      else
       if(iesup.gt.0) then  ! magnetization in the plane
        do ip=1,nh
         matrix(ip+nh,ip)=matrix(ip+nh,ip)-dup(1)*ss(ip)
         matrix(ip,ip+nh)=matrix(ip,ip+nh)-dup(1)*ss(ip)
        enddo
       endif
      endif

      if(iopt.ge.0) then 
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

      if(ieskin.ge.1) then
       if(iopt.ge.0) then 
        do ip=1,nh
         matrix(ip,ip)=matrix(ip,ip)-dek(1)
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)+dek(1)
        enddo
       else
        do ip=1,nh
         matrix(ip,ip)=matrix(ip,ip)-dek(1)
         matrix(ip+nh,ip+nh)=matrix(ip+nh,ip+nh)-dek(1)
        enddo
       endif
      endif

c Matrix diagonalization

      n1=4*nh*nh+1
      n2=n1+2*nh
      n3=n2+20*nh
      call dcopy(4*nh*nh,matrix,1,psip,1)

      call dsyev('V','L',2*nh,psip,2*nh,psip(n1),psip(n2),20*nh,info)

      do i=1,2*nh
       call dcopy(2*nh,psip(2*nh*(i-1)+1),1,w(1,i),1)
       eigenval(i)=psip(n1+i-1)
      enddo

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

      subroutine makeivic(nh,ixc,ivic)
      implicit none
      integer*4 i,nh
      integer*2 ivic(nh,*),ixc(nh,*)

      do i=1,nh
       if(ixc(i,2).eq.1) then
        ivic(i,2)=0
        ivic(i,3)=0
        ivic(i,4)=0
       endif
       if(ixc(i,2).eq.0) then
        ivic(i,6)=0
        ivic(i,7)=0
        ivic(i,8)=0
       endif
      enddo

      return
      end


