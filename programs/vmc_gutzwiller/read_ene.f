      program erread
      implicit none
      integer maxk,k,kt,i,j,nbin,npm,iskip,maxj,ibinit,icount
     >,ibin,nmis,lbin,nbuf
      parameter(nbuf=100,npm=10000)
      real*8 wsk,ebuf(nbuf),wbuf(nbuf),wtot(nbuf),wbufm(nbuf),eskip(npm)
     >,ek(0:nbuf),wk(0:nbuf),ebin(0:nbuf),wbin(0:nbuf),ebin2(0:nbuf)
            
      write(6,*) ' max k correction, bin lenght,bin init, iskip?' 
      read(5,*) maxk,lbin,ibinit,iskip
       
      open(unit=12,file='fort.12',form='unformatted',status='old')

      if(iskip.gt.npm-1) then
       write(6,*) 'recompile with larger npm',iskip
       stop
      endif
      if(iskip.eq.0) then 
       write(20,*) ' Calculation gs energy '
      endif

      open(unit=20,file='fort.20',form='formatted',status='unknown')

      rewind(12) 

      do i=0,nbuf
       ek(i)=0.d0
       wk(i)=0.d0
       ebin(i)=0.d0
       wbin(i)=0.d0
       ebin2(i)=0.d0
      enddo
      do i=1,nbuf
       wbuf(i)=0.d0
       wbufm(i)=0.d0
       wtot(i)=0.d0
       ebuf(i)=0.d0
      enddo

      icount=0
      ibin=0
      kt=0
      do while(kt.ge.0) 
       kt=kt+1
       do i=1,nbuf
        wbufm(i)=wbuf(i)
       enddo
       maxj=nbuf
       do j=1,maxj
        if(iskip.ne.0) then
         read(12,end=100) wbuf(j),wtot(j),(eskip(k),k=1,iskip),ebuf(j)
        else
         read(12,end=100) wbuf(j),wtot(j),ebuf(j)
        endif
       enddo

       if(kt.eq.1) then      
        icount=maxk  
        do j=maxk+1,maxj
         icount=icount+1
         do k=0,maxk
          if(k.eq.0)then
           wsk=wtot(j) 
          else
           wsk=wsk*wbuf(j-k)                     
          endif
          ek(k)=ek(k)+ebuf(j)*wsk
          wk(k)=wk(k)+wsk
         enddo   
         if(mod(icount,lbin).eq.0) then
          ibin=ibin+1
          if(ibin.ge.ibinit) then
           do k=0,maxk
            if(wk(k).ne.0.d0) then 
             ebin(k)=ebin(k)+ek(k)
             wbin(k)=wbin(k)+wk(k)
             ebin2(k)=ebin2(k)+ek(k)**2/wk(k)
            endif
           enddo
           if(wk(maxk).ne.0.d0) then 
            write(21,*) ek(maxk)/wk(maxk),wk(maxk),ibin 
           else
            write(21,*) 0.,wk(maxk),ibin 
           endif
          endif 
          do k=0,maxk
           ek(k)=0.d0
           wk(k)=0.d0
          enddo
         endif
        enddo   ! enddo j
       else
        do j=1,maxj
         icount=icount+1
         do k=0,maxk                              
          if(k.eq.0)then
           wsk=wtot(j)
          else
           if(j-k.ge.1) then
            wsk=wsk*wbuf(j-k)
           else
            wsk=wsk*wbufm(j-k+nbuf)                     
           endif
          endif 
          ek(k)=ek(k)+ebuf(j)*wsk
          wk(k)=wk(k)+wsk
         enddo
         if(mod(icount,lbin).eq.0) then
          ibin=ibin+1
          if(ibin.ge.ibinit) then
           do k=0,maxk
            if(wk(k).ne.0.d0) then 
             ebin(k)=ebin(k)+ek(k)
             wbin(k)=wbin(k)+wk(k)
             ebin2(k)=ebin2(k)+ek(k)**2/wk(k)
            endif
           enddo
           if(wk(maxk).ne.0.d0) then 
            write(21,*) ek(maxk)/wk(maxk),wk(maxk),ibin 
           else
            write(21,*) 0.,wk(maxk),ibin 
           endif
          endif 
          do k=0,maxk
           ek(k)=0.d0
           wk(k)=0.d0
          enddo
         endif 
        enddo   ! enddo j
       endif 
      enddo     ! enddo kt

100   continue
      nmis=ibin-ibinit+1
      write(6,*) ' number of measures done =',nmis
      do i=0,maxk
       ebin(i)=ebin(i)/wbin(i)
       ebin2(i)=dsqrt(dabs(ebin2(i)/wbin(i)-ebin(i)**2))
       ebin2(i)=ebin2(i)/dsqrt(dfloat(nmis))
      enddo
      write(20,*) ' Independent bins ',nmis,'of lenght ',lbin
      write(20,*) 
      write(20,*) ' Energy, error, # of bias correcting factors'
      do i=0,maxk
       write(20,*) ebin(i),ebin2(i),i
      enddo
      close(12)
      close(20) 
      stop
      end
