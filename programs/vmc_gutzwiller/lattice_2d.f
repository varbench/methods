      program lattice
      implicit none
      INTEGER(4) nbase,nsite,ntotalsite,nsimmax,nmom
      INTEGER(4) i,j,k,isg1,isg2,itest,ndim,n1,n2,iflag
     >,itest1,itest2,ii,it,jj,jt,l,m,nind,nindm,nindms
     >,maxmulti,maxmultis,ib1,ib2,iflagibc,ik1,ik2,is,ikm1,ikm2
     >,ikm,ik,mm,l1,l2,n1min,n1max,n2min,n2max,iflagref
     >,jshift,idists,ibroken,icount,jsav,mt
     >,iclass,index,iii,info,sumwar
      CHARACTER(5) oper1,oper2,oper3,oper4
      REAL(8) a1,a2,phi,pi,small,axx,ayy
     >,xi,yi,xj,yj,xij,yij,d2,xt,yt,d2t,det
     >,dij,d,test,test1,test2,uno,eig1,eig2,qx,qy
     >,phiKIN_r,phiKIN_i,phiBCS_r,phiBCS_i

      INTEGER(4), dimension(:), allocatable :: i1,i2
      INTEGER(4), dimension(:), allocatable :: imulti,imultis
      INTEGER(4), dimension(:), allocatable :: imultis_scra
      INTEGER(4), dimension(:), allocatable :: imark
      INTEGER(4), dimension(2) :: ipiv
      INTEGER(4), dimension(2,2) :: enne
      INTEGER(4), dimension(:), allocatable :: iwar
      INTEGER(4), dimension(:,:), allocatable :: warning
      INTEGER(4), dimension(:,:), allocatable :: boundary
      INTEGER(4), dimension(:,:), allocatable :: isymt
      INTEGER(4), dimension(:,:), allocatable :: imark2
      INTEGER(4), dimension(:,:,:), allocatable :: ivic,ivict,ivicns
      INTEGER(4), dimension(:,:,:), allocatable :: ivicns_scra
      REAL(8), dimension(:), allocatable :: xb,yb
      REAL(8), dimension(:), allocatable :: dindip
      REAL(8), dimension(2) :: Bbig
      REAL(8), dimension(2,2) :: Abig
      REAL(8), dimension(2,2) :: a,b,enne1,amat
      REAL(8), dimension(:,:), allocatable :: x,y
      REAL(8), dimension(:,:), allocatable :: dist
      COMPLEX(8), dimension(:,:,:), allocatable :: phaserKIN,phaserBCS
      COMPLEX(8), dimension(:,:,:), allocatable :: phaseiKIN,phaseiBCS
      COMPLEX(8), dimension(:,:,:), allocatable :: phaserKIN_scra
      COMPLEX(8), dimension(:,:,:), allocatable :: phaserBCS_scra
      COMPLEX(8), dimension(:,:,:), allocatable :: phaseiKIN_scra
      COMPLEX(8), dimension(:,:,:), allocatable :: phaseiBCS_scra

      open(unit=90,file='geometry.d',status='old')               !input file

      pi=dacos(-1.d0)
      small=1.d-7

c     ENTRIES      
c     a1=modulus of the first primitive vector
c     a2=modulus of the second primitive vector
c     phi=angle between the two primitive vectors
c
c     enne(i,j) is a 2x2 integer matrix defining the two translations which 
c     identify the cluster:
c     T_i = sum_j enne(i,j)*a_j   
c     We assume that the leading components are enne(1,1)>0 and enne(2,2)>0 

ccccccccccccccccccccccccc  READING PART ccccccccccccccccccccccccc
      read(90,*) oper1,a1,oper2,a2,phi
      if(oper1.eq.'sqrt') then
       a1=dsqrt(a1)
      elseif(oper1.ne.'') then
       write(6,*)'wrong operator'
       stop
      endif
      if(oper2.eq.'sqrt') then
       a2=dsqrt(a2)
      elseif(oper2.ne.'') then
       write(6,*)'wrong operator'
       stop
      endif
      phi=phi*pi/180.d0

      read(90,*) nbase                      !number of elements of the basis

      ALLOCATE(xb(nbase))
      ALLOCATE(yb(nbase))

      if(nbase.ne.1)then
       do i=1,nbase
        read(90,*) oper3,axx,oper4,ayy
        if(oper3.eq.'sqrt') axx=dsqrt(axx)
        if(oper4.eq.'sqrt') ayy=dsqrt(ayy)
        xb(i)=axx
        yb(i)=ayy
       enddo
      else
       xb(1)=0.d0
       yb(1)=0.d0
      endif

      read(90,*) enne(1,1),enne(1,2),enne(2,1),enne(2,2)
      read(90,*) ib1,ib2

      if(enne(1,1).le.0.or.enne(2,2).le.0)then
       write(6,*)'box positioning'
       stop
      endif

      nsite=abs(enne(1,1)*enne(2,2)-enne(1,2)*enne(2,1))

c     Primitive vectors
      a(1,1)=a1                !a_1 x component
      a(2,1)=0.d0              !a_1 y component
      a(1,2)=a2*dcos(phi)      !a_2 x component
      a(2,2)=a2*dsin(phi)      !a_2 y component
c     Reciprocal vectors
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      b(1,1)=a(2,2)/det
      b(2,1)=-a(1,2)/det
      b(1,2)=-a(2,1)/det
      b(2,2)=a(1,1)/det
ccccccccccccccccccccccccc  END READING PART ccccccccccccccccccccccccc

ccccccccccccccccccccccccc  CLUSTER PART ccccccccccccccccccccccccc
      ndim=0

      ALLOCATE(x(nsite,nbase))
      ALLOCATE(y(nsite,nbase))
      ALLOCATE(i1(nsite))
      ALLOCATE(i2(nsite))

      n1max=enne(1,1)+max(0,enne(2,1))
      n1min=min(0,enne(2,1))
      n2max=enne(2,2)+max(0,enne(1,2))
      n2min=min(0,enne(1,2))

      do i=1,2
       do j=1,2
        amat(i,j)=0.d0
        do k=1,2
         amat(i,j)=amat(i,j)+a(k,i)*a(k,j)       !a_i scalar a_j
        enddo
       enddo
      enddo

      do n1=n1min,n1max
       do n2=n2min,n2max
        do i=1,2
         do j=1,2
          Abig(i,j)=0.d0
          do k=1,2
           Abig(i,j)=Abig(i,j)+enne(j,k)*amat(k,i)
          enddo
         enddo
         Bbig(i)=n1*amat(1,i)+n2*amat(2,i)
        enddo
        call dgesv(2,1,Abig,2,ipiv,Bbig,2,info)
        if(info.ne.0) write(6,*) 'Info=',info
        iflag=0
        do i=1,2
         if(dabs(Bbig(i)).lt.small) Bbig(i)=0.d0
         if(dabs(Bbig(i)-1.d0).lt.small) Bbig(i)=1.d0
         if(Bbig(i).lt.0.d0.or.Bbig(i).ge.1.d0) iflag=1
        enddo
        if(iflag.eq.0)then
         ndim=ndim+1
         do j=1,nbase
          x(ndim,j)=xb(j)+dble(n1)*a(1,1)+dble(n2)*a(1,2)
          y(ndim,j)=yb(j)+dble(n1)*a(2,1)+dble(n2)*a(2,2)
         enddo
         i1(ndim)=n1
         i2(ndim)=n2
        endif
       enddo
      enddo

      if(nsite.ne.ndim)then
       write(6,*)'wrong site number',nsite,ndim
       stop
      endif

      write(6,*)'Nsite=',ndim
      write(6,*)'Nbase=',nbase

      write(6,*)
      write(6,*)'CLUSTER'
      write(6,*)
      do i=1,ndim
       do j=1,nbase
        write(6,*) (i-1)*nbase+j,x(i,j),y(i,j)
        write(88,*) (i-1)*nbase+j,x(i,j),y(i,j)
       enddo
      enddo

c  Distances

      ntotalsite=nsite*nbase

      ALLOCATE(dist(ntotalsite,ntotalsite))
      ALLOCATE(boundary(ntotalsite,ntotalsite))
      ALLOCATE(warning(ntotalsite,ntotalsite))

c     write(6,*)
c     write(6,*)'DISTANCES'
c     write(6,*)

      do i=1,ndim
       do ii=1,nbase
        it=(i-1)*nbase+ii
        xi=x(i,ii)
        yi=y(i,ii)
        do j=1,ndim
         do jj=1,nbase
          jt=(j-1)*nbase+jj
          xj=x(j,jj)
          yj=y(j,jj)
          xij=xi-xj
          yij=yi-yj
          d2=xij**2+yij**2
          dist(it,jt)=d2
          iflagibc=0
          do isg1=-1,1
           do isg2=-1,1
            xt=xij+
     &      (isg1*enne(1,1)+isg2*enne(2,1))*a(1,1)+
     &      (isg1*enne(1,2)+isg2*enne(2,2))*a(1,2)
            yt=yij+
     &      (isg1*enne(1,1)+isg2*enne(2,1))*a(2,1)+
     &      (isg1*enne(1,2)+isg2*enne(2,2))*a(2,2)
            d2t=xt**2+yt**2
            if(d2t.lt.dist(it,jt)) then
             dist(it,jt)=d2t
             iflagibc=0
             if(mod(ib1*abs(isg1)+ib2*abs(isg2),2).ne.0)
     >       iflagibc=1
            endif
           enddo
          enddo
          dist(it,jt)=dsqrt(dist(it,jt))
          boundary(it,jt)=1
          if(iflagibc.eq.1) boundary(it,jt)=-1
c         write(6,*) it,jt,dist(it,jt)
         enddo
        enddo
       enddo
      enddo

      do it=1,ntotalsite
       do jt=1,ntotalsite
        warning(it,jt)=0
       enddo
      enddo

      do i=1,ndim
       do ii=1,nbase
        it=(i-1)*nbase+ii
        xi=x(i,ii)
        yi=y(i,ii)
        do j=1,ndim
         do jj=1,nbase
          jt=(j-1)*nbase+jj
          xj=x(j,jj)
          yj=y(j,jj)
          xij=xi-xj
          yij=yi-yj
          iflagibc=0
          do isg1=-1,1
           do isg2=-1,1
            xt=xij+
     &      (isg1*enne(1,1)+isg2*enne(2,1))*a(1,1)+
     &      (isg1*enne(1,2)+isg2*enne(2,2))*a(1,2)
            yt=yij+
     &      (isg1*enne(1,1)+isg2*enne(2,1))*a(2,1)+
     &      (isg1*enne(1,2)+isg2*enne(2,2))*a(2,2)
            d2t=dsqrt(xt**2+yt**2)
            if(dabs(d2t-dist(it,jt)).lt.small) then
             iflagibc=0
             if(mod(ib1*abs(isg1)+ib2*abs(isg2),2).ne.0)
     >       iflagibc=1
             if((iflagibc.eq.0.and.boundary(it,jt).eq.-1).or.
     >          (iflagibc.eq.1.and.boundary(it,jt).eq.1)) then
              warning(it,jt)=1
             endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

c  independent distances

      ALLOCATE(dindip(ntotalsite))

      do i=1,nsite*nbase
       dindip(i)=0.d0
      enddo

      dindip(1)=dist(1,1)
      k=1
      do ii=1,nbase
       do j=ii+1,ndim*nbase
        do l=1,k
         test=dist(ii,j)-dindip(l)
         if(abs(test).lt.small) goto 21
         if(test.lt.0.d0) goto 11
        enddo
        l=k+1
11      continue            ! j is at l-th position
        if(l.le.k)then
         do m=k,l,-1
          dindip(m+1)=dindip(m)
         enddo
        endif
        dindip(l)=dist(ii,j)
        k=k+1
21      continue            ! j is already in the list
       enddo
      enddo

      do i=2,k
       dindip(i-1)=dindip(i)
      enddo
      dindip(k)=0.d0

      nind=k-1               ! number of independent distances

      write(6,*)
      write(6,*) 'TOTAL INDEPENDENT DISTANCES   =',nind
      write(6,*)
      do i=1,nind
       write(6,*) i,dindip(i)
      enddo

c  inverse N matrix
      det=enne(1,1)*enne(2,2)-enne(1,2)*enne(2,1)
      enne1(1,1)=enne(2,2)/det
      enne1(1,2)=-enne(1,2)/det
      enne1(2,1)=-enne(2,1)/det
      enne1(2,2)=enne(1,1)/det

c  Allowed momenta (in units of 2*pi)

      nmom=0
      n1max=enne(1,1)+max(0,enne(1,2))
      n1min=min(0,enne(1,2))
      n2max=enne(2,2)+max(0,enne(2,1))
      n2min=min(0,enne(2,1))

      uno=1.d0-small
      do n1=n1min,n1max
       do n2=n2min,n2max
        eig1=enne1(1,1)*n1+enne1(1,2)*n2
        eig2=enne1(2,1)*n1+enne1(2,2)*n2
        if(eig1.ge.0.d0.and.eig1.lt.uno.and.
     >     eig2.ge.0.d0.and.eig2.lt.uno)then
         qx=eig1*b(1,1)+eig2*b(1,2)
         qy=eig1*b(2,1)+eig2*b(2,2)
         nmom=nmom+1
         write(27,*) qx,qy
        endif
       enddo
      enddo

      write(6,*)'Independent momenta=',nmom

      if(nmom.ne.ndim)then
       write(6,*)'problems in determining the momenta'
       stop
      endif

      do i=1,ndim
       write(28,*) x(i,1),y(i,1)
      enddo

ccccccccccccccccccccccccc  END CLUSTER PART ccccccccccccccccccccccccc

ccccccccccccccccccccccccc  TABLE OF NEIGHBORS ccccccccccccccccccccccccc

      nindm=nind
      ALLOCATE(imulti(nindm))
      ALLOCATE(imark(ndim*nbase))

      do k=1,nindm
       d=dindip(k)
       do i=1,ndim*nbase
        imark(i)=0
        do j=1,ndim*nbase
         dij=dist(i,j)
         if(dabs(dij-d).lt.small) imark(i)=imark(i)+1
        enddo
       enddo
       maxmulti=0
       do i=1,ndim*nbase
        maxmulti=max(maxmulti,imark(i))
       enddo
       imulti(k)=maxmulti
      enddo

      DEALLOCATE(imark)

      maxmulti=0
      do k=1,nindm
       maxmulti=max(maxmulti,imulti(k))
      enddo

      ALLOCATE(ivic(ntotalsite,maxmulti,nindm))
      ALLOCATE(iwar(nindm))

      do k=1,nindm
       iwar(k)=0
       do i=1,ndim*nbase
        do j=1,maxmulti
         ivic(i,j,k)=0
        enddo
       enddo
      enddo

      do k=1,nindm
       d=dindip(k)
       do i=1,ndim*nbase
        ii=0
        do j=1,ndim*nbase
         dij=dist(i,j)
         if(abs(dij-d).lt.small) then
          ii=ii+1
          ivic(i,ii,k)=j
          if(boundary(i,j).lt.0) ivic(i,ii,k)=-ivic(i,ii,k)
          if(warning(i,j).eq.1) iwar(k)=1
         endif
        enddo
       enddo
      enddo

c  In this table the order of neighbors for each site is given by the 
c  lexicographic order of the sites
c     write(6,*) 
c     write(6,*) 'Table of neighbors'
c     do k=1,nindm
c      write(6,*) 'At distance k=',k
c      do i=1,ndim*nbase
c       write(6,334) i,(ivic(i,j,k),j=1,imulti(k))
c      enddo
c     enddo

      ALLOCATE(isymt(ntotalsite,ndim))  ! table of translational symmetries

      do i=1,ndim            ! sites
       it=0
       do k=1,ndim           ! translations
        ik1=i1(i)+i1(k)
        ik2=i2(i)+i2(k)
        it=it+1
        iflag=0
        do l=1,ndim
         l1=ik1-i1(l)
         l2=ik2-i2(l)
         test1=enne1(1,1)*l1+enne1(2,1)*l2
         test2=enne1(1,2)*l1+enne1(2,2)*l2
         if(itest(test1).eq.0.and.itest(test2).eq.0) then
          iflag=iflag+1
          ik=l
         endif
        enddo
        if(iflag.ne.1) then
         write(6,*)'wrong reduction'
         stop
        endif
        do ii=1,nbase
         j=(i-1)*nbase+ii
         jj=(ik-1)*nbase+ii
         isymt(j,it)=jj
        enddo
       enddo
      enddo

      ALLOCATE(ivict(ntotalsite,maxmulti,nindm))
      ALLOCATE(imark(ndim))

      do k=1,nindm
       do i=1,ndim*nbase
        do j=1,maxmulti
         ivict(i,j,k)=0
        enddo
       enddo
      enddo

      do i=1,ndim
       imark(i)=0
      enddo

      IF(nbase.ne.1) THEN ! we do not use inversion symmetry
       do k=1,nindm
        do i=1,nbase
         do j=1,imulti(k) 
          ii=abs(ivic(i,j,k))
          if(ii.ne.0) then
           do it=1,nsite   ! traslations
            m=isymt(i,it)
            l=isymt(ii,it)
            ivict(m,j,k)=l
            if(boundary(m,l).lt.0) ivict(m,j,k)=-ivict(m,j,k)
           enddo
          endif
         enddo
        enddo
       enddo
      ELSE                ! we order neighbors using the inversion symmetry 
       do k=1,nindm
        do ii=1,nsite
         imark(ii)=0
        enddo
        jshift=imulti(k)/2
        jj=0
        do j=1,imulti(k) 
         i=abs(ivic(1,j,k))
         if(imark(i).eq.0) then
          imark(i)=1  ! mark the site...
          iflag=0
          ik1=i1(i)-i1(1)
          ik2=i2(i)-i2(1)
          ikm1=-ik1+i1(1)
          ikm2=-ik2+i2(1)
          do l=1,ndim              !reduction
           l1=ikm1-i1(l)
           l2=ikm2-i2(l)
           test1=enne1(1,1)*l1+enne1(2,1)*l2
           test2=enne1(1,2)*l1+enne1(2,2)*l2
           if(itest(test1).eq.0.and.itest(test2).eq.0) then
            iflag=iflag+1
            ikm=l         ! inverted site
           endif
          enddo
          if(iflag.ne.1) then
           write(6,*) 'wrong reduction'
           stop
          endif
          if(i.ne.ikm) then ! there is an inverted site
           jj=jj+1
           iflagref=0
           imark(ikm)=1  ! ...and its inverted
          elseif(i.eq.ikm) then ! the inverted site is the site itself
           jj=jj+1
           iflagref=1
          endif
          do it=1,nsite   ! traslations
           m=isymt(1,it)
           l=isymt(i,it)
           if(iflagref.eq.0) then
            ivict(m,jj,k)=l
            if(boundary(m,l).lt.0) ivict(m,jj,k)=-ivict(m,jj,k)
            l=isymt(ikm,it)
            ivict(m,jj+jshift,k)=l
            if(boundary(m,l).lt.0) 
     >      ivict(m,jj+jshift,k)=-ivict(m,jj+jshift,k)
           elseif(iflagref.eq.1) then
            ivict(m,jj,k)=l
            if(boundary(m,l).lt.0)
     >      ivict(m,jj,k)=-ivict(m,jj,k)
           endif
          enddo
         endif
        enddo
       enddo
      ENDIF

c  In this table the order of neighbors for each site is given by the 
c  order of the first site
      write(6,*) 
      write(6,*) 'Ordered Table of neighbors'
      do k=1,nindm
       write(6,*) 'At distance k=',k
       do i=1,ndim*nbase
        write(6,334) i,(ivict(i,j,k),j=1,imulti(k))
       enddo
      enddo

334   format(i4,4x,100(i5,2x))

c  Table with phases

      write(6,*) 'How many classes for the WF?'
      read(5,*) idists
      if(idists.gt.nindm) then
       write(6,*) 'Too many classes',idists,nindm
       stop
      endif
      sumwar=0
      do i=1,idists
       sumwar=sumwar+iwar(i)
      enddo
      if(sumwar.gt.0) then
       write(6,*) 'Distances not compatible with boundary conditions'
       stop
      endif
      write(6,*) 'Do you want to break symmetries inside the unit cell?'
      write(6,*) 'NO = 0, YES = 1'
      read(5,*) ibroken

      nindms=idists
      if(ibroken.eq.1) then
       icount=0
       do k=1,idists
        ALLOCATE(imark2(nbase,imulti(k)))
        do ii=1,nbase
         do j=1,imulti(k)
          imark2(ii,j)=0
         enddo
        enddo
        do j=1,imulti(k)
         do ii=1,nbase
          if(imark2(ii,j).eq.0) then
           imark2(ii,j)=1
           i=abs(ivict(ii,j,k))
           if(i.ne.0) then
            itest1=(i-1)/nbase+1
            itest2=i-(itest1-1)*nbase
            if(nbase.eq.1) then
             iflagref=0                    ! there is the inverted
             if(mod(imulti(k),2).ne.0.and.j.eq.imulti(k)) 
     >       iflagref=1                    ! there is no inverted
             if(iflagref.eq.0) then
              jshift=imulti(k)/2
              imark2(ii,j+jshift)=1        ! mark the inverted bond
             endif
            elseif(nbase.gt.1) then
             jsav=0
             do jj=1,imulti(k)
              if(abs(ivict(i,jj,k)).eq.ii) jsav=jj
             enddo
             if(jsav.eq.0) then
              write(6,*) 'Something went wrong!'
              stop
             endif
             imark2(itest2,jsav)=1 ! mark the opposite bond
            endif 
            icount=icount+1
           endif
          endif
         enddo
        enddo
        DEALLOCATE(imark2)
       enddo
       nindms=icount
      endif

      write(6,*)
      write(6,*) 'Number of parameters=',nindms
      write(6,*)

      IF(ibroken.eq.0) THEN
       maxmultis=maxmulti
       ALLOCATE(imultis(nindms))
       ALLOCATE(ivicns(ntotalsite,maxmultis,nindms))
       ALLOCATE(phaserKIN(ntotalsite,maxmultis,nindms))
       ALLOCATE(phaseiKIN(ntotalsite,maxmultis,nindms))
       ALLOCATE(phaserBCS(ntotalsite,maxmultis,nindms))
       ALLOCATE(phaseiBCS(ntotalsite,maxmultis,nindms))
       do k=1,nindms
        imultis(k)=0
        do j=1,maxmultis
         do ii=1,nsite*nbase
          ivicns(ii,j,k)=0
          phaserKIN(ii,j,k)=dcmplx(0.d0,0.d0)
          phaseiKIN(ii,j,k)=dcmplx(0.d0,0.d0)
          phaserBCS(ii,j,k)=dcmplx(0.d0,0.d0)
          phaseiBCS(ii,j,k)=dcmplx(0.d0,0.d0)
         enddo
        enddo
       enddo
       do k=1,idists
        imultis(k)=imulti(k)
        do j=1,imultis(k)
         do ii=1,nsite*nbase
          ivicns(ii,j,k)=ivict(ii,j,k)
         enddo
        enddo
        ALLOCATE(imark2(nbase,imultis(k)))
        do ii=1,nbase
         do j=1,imultis(k)
          imark2(ii,j)=0
         enddo
        enddo
        do ii=1,nbase
         do j=1,imultis(k)
          if(imark2(ii,j).eq.0) then
           imark2(ii,j)=1
           i=abs(ivict(ii,j,k))
           if(i.ne.0) then
            itest1=(i-1)/nbase+1
            itest2=i-(itest1-1)*nbase
            if(nbase.eq.1) then
             iflagref=0                     ! there is the inverted
             if(mod(imultis(k),2).ne.0.and.j.eq.imultis(k)) 
     >       iflagref=1                     ! there is no inverted
             if(iflagref.eq.0) then
              jshift=imultis(k)/2
              imark2(ii,j+jshift)=1         ! mark the inverted bond
             endif
            elseif(nbase.gt.1) then
             jsav=0
             do jj=1,imultis(k)
              if(abs(ivict(i,jj,k)).eq.ii) jsav=jj
             enddo
             if(jsav.eq.0) then
              write(6,*) 'Something went wrong!'
              stop
             endif
             imark2(itest2,jsav)=1 ! mark the opposite bond
            endif 
            write(6,*) ' PARAMETER k=',k
            write(6,*) ' Choose phases for the bond',ii,i
            write(6,*) ' hopping (real and imaginary) and pairing
     > (real and imaginary) in degrees'
            read(5,*) phiKIN_r,phiKIN_i,phiBCS_r,phiBCS_i
            phiKIN_r=pi*phiKIN_r/180.d0
            phiKIN_i=pi*phiKIN_i/180.d0
            phiBCS_r=pi*phiBCS_r/180.d0
            phiBCS_i=pi*phiBCS_i/180.d0
            do it=1,nsite   ! traslations
             m=isymt(ii,it)
             phaserKIN(m,j,k)=dcmplx(dcos(phiKIN_r),dsin(phiKIN_r))
             phaseiKIN(m,j,k)=dcmplx(dcos(phiKIN_i),dsin(phiKIN_i))
             phaserBCS(m,j,k)=dcmplx(dcos(phiBCS_r),dsin(phiBCS_r))
             phaseiBCS(m,j,k)=dcmplx(dcos(phiBCS_i),dsin(phiBCS_i))
             if(nbase.eq.1.and.iflagref.eq.0) then
              phaserKIN(m,j+jshift,k)=dconjg(phaserKIN(m,j,k))
              phaseiKIN(m,j+jshift,k)=-dconjg(phaseiKIN(m,j,k))
              phaserBCS(m,j+jshift,k)=phaserBCS(m,j,k)
              phaseiBCS(m,j+jshift,k)=phaseiBCS(m,j,k)
             elseif(nbase.gt.1) then
              mt=isymt(itest2,it)
              phaserKIN(mt,jsav,k)=dconjg(phaserKIN(m,j,k))
              phaseiKIN(mt,jsav,k)=-dconjg(phaseiKIN(m,j,k))
              phaserBCS(mt,jsav,k)=phaserBCS(m,j,k)
              phaseiBCS(mt,jsav,k)=phaseiBCS(m,j,k)
             endif
            enddo
           endif
          endif
         enddo
        enddo
        DEALLOCATE(imark2)
       enddo
      ELSEIF(ibroken.eq.1) THEN
c  Here you can put the constraint that different bonds belong to the same class
c  All different bonds are listed and you have the option to choose if
c  a bond defines a new class or it is equal to a previuosly defined class.
c  see the variable ICLASS
       ALLOCATE(imultis_scra(nindms))
       ALLOCATE(ivicns_scra(ntotalsite,100,nindms))
       ALLOCATE(phaserKIN_scra(ntotalsite,100,nindms))
       ALLOCATE(phaseiKIN_scra(ntotalsite,100,nindms))
       ALLOCATE(phaserBCS_scra(ntotalsite,100,nindms))
       ALLOCATE(phaseiBCS_scra(ntotalsite,100,nindms))
       do k=1,nindms
        imultis_scra(k)=0
        do j=1,100
         do ii=1,nsite*nbase
          ivicns_scra(ii,j,k)=0
          phaserKIN_scra(ii,j,k)=dcmplx(0.d0,0.d0)
          phaseiKIN_scra(ii,j,k)=dcmplx(0.d0,0.d0)
          phaserBCS_scra(ii,j,k)=dcmplx(0.d0,0.d0)
          phaseiBCS_scra(ii,j,k)=dcmplx(0.d0,0.d0)
         enddo
        enddo
       enddo
       icount=0
       do k=1,idists
        ALLOCATE(imark2(nbase,imulti(k)))
        do ii=1,nbase
         do j=1,imulti(k)
          imark2(ii,j)=0
         enddo
        enddo
        do j=1,imulti(k)
         do ii=1,nbase
          if(imark2(ii,j).eq.0) then
           imark2(ii,j)=1
           i=abs(ivict(ii,j,k))
           if(i.ne.0) then
            itest1=(i-1)/nbase+1
            itest2=i-(itest1-1)*nbase
            if(nbase.eq.1) then
             iflagref=0                    ! there is the inverted
             if(mod(imulti(k),2).ne.0.and.j.eq.imulti(k)) 
     >       iflagref=1                    ! there is no inverted
             if(iflagref.eq.0) then
              jshift=imulti(k)/2
              imark2(ii,j+jshift)=1        ! mark the inverted bond
             endif
            elseif(nbase.gt.1) then
             jsav=0
             do jj=1,imulti(k)
              if(abs(ivict(i,jj,k)).eq.ii) jsav=jj
             enddo
             if(jsav.eq.0) then
              write(6,*) 'Something went wrong!'
              stop
             endif
             imark2(itest2,jsav)=1 ! mark the opposite bond
            endif 
            write(6,*) 'Sites=',ii,i
            write(6,*) 'Do you want to consider it as a NEW bond?'
            write(6,*) 'Yes: put 0'
            write(6,*) 'NO : put the number of the class <=',icount 
111         continue
            read(5,*) iclass
            if(iclass.gt.icount) then
             write(6,*) 'Class not yet defined!'
             goto 111
            endif 
            if(iclass.eq.0) then
             icount=icount+1
             index=icount
             write(6,*) 
             write(6,335) index,ii,i
             write(6,*) 
            else
             index=iclass
             write(6,*) 
             write(6,335) index,ii,i
             write(6,*) 
            endif
            imultis_scra(index)=imultis_scra(index)+1
            iii=imultis_scra(index)
            phiKIN_r=0.d0
            phiKIN_i=0.d0
            phiBCS_r=0.d0
            phiBCS_i=0.d0
c            if(iclass.ne.0) then
             write(6,*) ' Choose phases for the bond',ii,i
             write(6,*) ' for hopping (real and imaginary) and pairing 
     >       (real and imaginary) in degrees'
             read(5,*) phiKIN_r,phiKIN_i,phiBCS_r,phiBCS_i
             phiKIN_r=pi*phiKIN_r/180.d0
             phiKIN_i=pi*phiKIN_i/180.d0
             phiBCS_r=pi*phiBCS_r/180.d0
             phiBCS_i=pi*phiBCS_i/180.d0
c            endif
            do it=1,nsite   ! traslations
             m=isymt(ii,it)
             ivicns_scra(m,iii,index)=ivict(m,j,k)
             phaserKIN_scra(m,iii,index)=
     >       dcmplx(dcos(phiKIN_r),dsin(phiKIN_r))
             phaseiKIN_scra(m,iii,index)=
     >       dcmplx(dcos(phiKIN_i),dsin(phiKIN_i))
             phaserBCS_scra(m,iii,index)=
     >       dcmplx(dcos(phiBCS_r),dsin(phiBCS_r))
             phaseiBCS_scra(m,iii,index)=
     >       dcmplx(dcos(phiBCS_i),dsin(phiBCS_i))
            enddo
            if(nbase.gt.1) then
             if(ii.ne.itest2) then
              do it=1,nsite   ! traslations
               mt=isymt(itest2,it)
               ivicns_scra(mt,iii,index)=ivict(mt,jsav,k)
               phaserKIN_scra(mt,iii,index)=
     >         dconjg(phaserKIN_scra(m,iii,index))
               phaseiKIN_scra(mt,iii,index)=
     >         -dconjg(phaseiKIN_scra(m,iii,index))
               phaserBCS_scra(mt,iii,index)=phaserBCS_scra(m,iii,index)
               phaseiBCS_scra(mt,iii,index)=phaseiBCS_scra(m,iii,index)
              enddo
             elseif(ii.eq.itest2) then
              imultis_scra(index)=imultis_scra(index)+1
              iii=imultis_scra(index)
              do it=1,nsite   ! traslations
               mt=isymt(ii,it)
               ivicns_scra(mt,iii,index)=ivict(mt,jsav,k)
               phaserKIN_scra(mt,iii,index)=
     >         dcmplx(dcos(phiKIN_r),-dsin(phiKIN_r))
               phaseiKIN_scra(mt,iii,index)=
     >         -dcmplx(dcos(phiKIN_i),-dsin(phiKIN_i))
               phaserBCS_scra(mt,iii,index)=
     >         dcmplx(dcos(phiBCS_r),dsin(phiBCS_r))
               phaseiBCS_scra(mt,iii,index)=
     >         dcmplx(dcos(phiBCS_i),dsin(phiBCS_i))
              enddo
             endif
            endif
            if(nbase.eq.1.and.iflagref.eq.0) then
             imultis_scra(index)=imultis_scra(index)+1
             iii=imultis_scra(index)
             do it=1,nsite   ! traslations
              m=isymt(ii,it)
              ivicns_scra(m,iii,index)=ivict(m,j+jshift,k)
              phaserKIN_scra(m,iii,index)=
     >        dcmplx(dcos(phiKIN_r),-dsin(phiKIN_r))
              phaseiKIN_scra(m,iii,index)=
     >        -dcmplx(dcos(phiKIN_i),-dsin(phiKIN_i))
              phaserBCS_scra(m,iii,index)=
     >        dcmplx(dcos(phiBCS_r),dsin(phiBCS_r))
              phaseiBCS_scra(m,iii,index)=
     >        dcmplx(dcos(phiBCS_i),dsin(phiBCS_i))
             enddo
            endif
           endif
          endif
         enddo
        enddo
        DEALLOCATE(imark2)
       enddo
       maxmultis=0
       do k=1,nindms
        maxmultis=max(maxmultis,imultis_scra(k))
       enddo
       nindms=icount
       ALLOCATE(imultis(nindms))
       ALLOCATE(ivicns(ntotalsite,maxmultis,nindms))
       ALLOCATE(phaserKIN(ntotalsite,maxmultis,nindms))
       ALLOCATE(phaseiKIN(ntotalsite,maxmultis,nindms))
       ALLOCATE(phaserBCS(ntotalsite,maxmultis,nindms))
       ALLOCATE(phaseiBCS(ntotalsite,maxmultis,nindms))
       do k=1,nindms
        imultis(k)=imultis_scra(k)
        do j=1,maxmultis
         do ii=1,nsite*nbase
          ivicns(ii,j,k)=ivicns_scra(ii,j,k)
          phaserKIN(ii,j,k)=phaserKIN_scra(ii,j,k)
          phaseiKIN(ii,j,k)=phaseiKIN_scra(ii,j,k)
          phaserBCS(ii,j,k)=phaserBCS_scra(ii,j,k)
          phaseiBCS(ii,j,k)=phaseiBCS_scra(ii,j,k)
         enddo
        enddo
       enddo
       DEALLOCATE(imultis_scra)
       DEALLOCATE(ivicns_scra)
       DEALLOCATE(phaserKIN_scra)
       DEALLOCATE(phaseiKIN_scra)
       DEALLOCATE(phaserBCS_scra)
       DEALLOCATE(phaseiBCS_scra)
      ENDIF

335   format(1x,'Class=',i3,4x,'site_1=',i4,3x,'site_2=',i4)

      if(idists.ne.0) then
       write(6,*) 
       write(6,*) 'Table of neighbors and phases'
       do k=1,nindms
        write(6,*) 'For the class k=',k
        do i=1,ndim*nbase
         write(6,336) i,(ivicns(i,j,k),j=1,imultis(k))
         write(6,337) i,(phaserKIN(i,j,k),j=1,imultis(k))
         write(6,337) i,(phaseiKIN(i,j,k),j=1,imultis(k))
         write(6,337) i,(phaserBCS(i,j,k),j=1,imultis(k))
         write(6,337) i,(phaseiBCS(i,j,k),j=1,imultis(k))
        enddo
       enddo
      endif

336   format(i4,3x,100(i5,7x))
337   format(i4,4x,100(f10.6,2x))

ccccccccccccccccccccccccc  END TABLE OF NEIGHBORS ccccccccccccccccccccccccc

ccccccccccccccccccccccccc  WRITE INFORMATION ccccccccccccccccccccccccc

      rewind(9)
      write(9) ndim,nbase,nindm,nindms,maxmulti,maxmultis
      write(9) imulti,imultis,ivict,ivicns
      write(9) phaserKIN,phaserBCS,phaseiKIN,phaseiBCS

      rewind(8)
      write(8) ndim,nbase
      write(8) isymt

ccccccccccccccccccccccccc  END WRITE INFORMATION ccccccccccccccccccccccccc

      DEALLOCATE(xb)
      DEALLOCATE(yb)
      DEALLOCATE(x)
      DEALLOCATE(y)
      DEALLOCATE(i1)
      DEALLOCATE(i2)
      DEALLOCATE(dist)
      DEALLOCATE(iwar)
      DEALLOCATE(warning)
      DEALLOCATE(boundary)
      DEALLOCATE(dindip)
      DEALLOCATE(ivic)
      DEALLOCATE(ivict)
      DEALLOCATE(imulti)
      DEALLOCATE(imultis)
      DEALLOCATE(isymt)
      DEALLOCATE(imark)
      DEALLOCATE(ivicns)
      DEALLOCATE(phaserKIN)
      DEALLOCATE(phaseiKIN)
      DEALLOCATE(phaserBCS)
      DEALLOCATE(phaseiBCS)

      stop
      end

      integer*4 function itest(test)
      implicit none
      integer*4 ii
      real*8 test,check,small
      parameter(small=1.d-10)
      
      ii=nint(test)
      check=test-dble(ii)
      if(abs(check).lt.small)then
       itest=0      
      else
       itest=1
      endif

      return
      end
