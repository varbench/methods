      program lattice
      implicit none
      INTEGER(4) nbase,nsite,ntotalsite,nsimmax,nmom,enne
      INTEGER(4) i,j,k,isg1,isg2,isg3,itest,ndim,n1,n2,n3,iflag
     >,itest1,itest2,ii,it,jj,jt,l,m,nind,nindm,nindms
     >,maxmulti,maxmultis,ib1,iflagibc,ik1,is
     >,ikm1,ikm2,ikm3,ikm,ik,mm,l1
     >,iflagref,jshift,idists,ibroken,icount,jsav,mt
     >,iclass,index,iii,info,sumwar
      CHARACTER(5) oper1,oper2
      REAL(8) a11,pi,small,axx,xi,xij,d2,xt,d2t,xj
     >,det,dij,d,test,test1,test2,test3,uno,qx,enne1
     >,phiKIN_r,phiKIN_i,phiBCS_r,phiBCS_i

      INTEGER(4), dimension(:), allocatable :: i1
      INTEGER(4), dimension(:), allocatable :: imulti,imultis
      INTEGER(4), dimension(:), allocatable :: imultis_scra
      INTEGER(4), dimension(:), allocatable :: imark
      INTEGER(4), dimension(:), allocatable :: iwar
      INTEGER(4), dimension(:,:), allocatable :: warning
      INTEGER(4), dimension(:,:), allocatable :: boundary
      INTEGER(4), dimension(:,:), allocatable :: isymt
      INTEGER(4), dimension(:,:), allocatable :: imark2
      INTEGER(4), dimension(:,:,:), allocatable :: ivic,ivict,ivicns
      INTEGER(4), dimension(:,:,:), allocatable :: ivicns_scra
      REAL(8), dimension(:), allocatable :: xb
      REAL(8), dimension(:), allocatable :: dindip
      REAL(8), dimension(:,:), allocatable :: x
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
c     a11=the length of the primitive vector
c
c     enne is the integer defining the translation which identifies the
cluster:

ccccccccccccccccccccccccc  READING PART ccccccccccccccccccccccccc
      read(90,*) oper1,a11
      if(oper1.eq.'sqrt') then
       a11=dsqrt(a11)
      elseif(oper1.ne.'') then
       write(6,*)'wrong operator'
       stop
      endif

      read(90,*) nbase                      !number of elements of the basis

      ALLOCATE(xb(nbase))

      if(nbase.ne.1)then
       do i=1,nbase
        read(90,*) oper2,axx
        if(oper2.eq.'sqrt') axx=dsqrt(axx)
        xb(i)=axx
       enddo
      else
       xb(1)=0.d0
      endif

      read(90,*) enne
      read(90,*) ib1

      if(enne.lt.0)then
       write(6,*)'box positioning'
       stop
      endif

      nsite=enne

ccccccccccccccccccccccccc  END READING PART ccccccccccccccccccccccccc

ccccccccccccccccccccccccc  CLUSTER PART ccccccccccccccccccccccccc
      ndim=0

      ALLOCATE(x(nsite,nbase))
      ALLOCATE(i1(nsite))

      ndim=0
      do n1=1,nsite
       ndim=ndim+1
       do j=1,nbase
        x(ndim,j)=xb(j)+a11*dble(n1-1)
       enddo
       i1(ndim)=n1-1
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
        write(6,666) (i-1)*nbase+j,x(i,j)
       enddo
      enddo
666   format(i4,f10.6)

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
        do j=1,ndim
         do jj=1,nbase
          jt=(j-1)*nbase+jj
          xj=x(j,jj)
          xij=xi-xj
          d2=xij**2
          dist(it,jt)=d2
          iflagibc=0
          do isg1=-1,1
           xt=xij+isg1*enne*a11
           d2t=xt**2
           if(d2t.lt.dist(it,jt)) then
            dist(it,jt)=d2t
            iflagibc=0
            if(mod(ib1*abs(isg1),2).ne.0) iflagibc=1
           endif
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
        do j=1,ndim
         do jj=1,nbase
          jt=(j-1)*nbase+jj
          xj=x(j,jj)
          xij=xi-xj
          iflagibc=0
          do isg1=-1,1
           xt=xij+isg1*enne*a11
           d2t=dsqrt(xt**2)
           if(dabs(d2t-dist(it,jt)).lt.small) then
            iflagibc=0
            if(mod(ib1*abs(isg1),2).ne.0) iflagibc=1
            if((iflagibc.eq.0.and.boundary(it,jt).eq.-1).or.
     >         (iflagibc.eq.1.and.boundary(it,jt).eq.1)) then
             warning(it,jt)=1
            endif
           endif
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
      enne1=1.d0/dble(enne)

c  Allowed momenta (in units of 2*pi)

      nmom=0
      do n1=0,enne-1
       nmom=nmom+1
       qx=enne1*n1
       write(27,*) qx
      enddo

      write(6,*)'Independent momenta=',nmom

      if(nmom.ne.ndim)then
       write(6,*)'problems in determining the momenta'
       stop
      endif

      do i=1,ndim
       write(28,*) x(i,1)
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
        it=it+1
        iflag=0
        do l=1,ndim
         l1=ik1-i1(l)
         test1=enne1*l1
         if(itest(test1).eq.0)then
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
          ikm1=-ik1+i1(1)
          do l=1,ndim              !reduction
           l1=ikm1-i1(l)
           test1=enne1*l1
           if(itest(test1).eq.0) then
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
            if(iclass.ne.0) then
            write(6,*) ' Choose phases for the bond',ii,i
             write(6,*) ' for hopping (real and imaginary) and pairing
     >       (real and imaginary) in degrees'
             read(5,*) phiKIN_r,phiKIN_i,phiBCS_r,phiBCS_i
             phiKIN_r=pi*phiKIN_r/180.d0
             phiKIN_i=pi*phiKIN_i/180.d0
             phiBCS_r=pi*phiBCS_r/180.d0
             phiBCS_i=pi*phiBCS_i/180.d0
            endif
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
      DEALLOCATE(x)
      DEALLOCATE(i1)
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
