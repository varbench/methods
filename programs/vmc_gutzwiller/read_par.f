      program readall
      implicit none
      integer*4 ngen,ibinit,np,nmat,npm,ibin,nmis,kk,kt
      parameter(npm=10000)
      real*8 ebuf(npm),ebin(npm),ebin2(npm),wbin

      open(unit=12,file='fort.13',form='unformatted',status='old')
      open(unit=20,file='fort.20',form='formatted',status='unknown')
      open(unit=21,file='fort.21',form='formatted',status='unknown')

      write(6,*) 'ibinit, np?'
      read(5,*) ibinit,np

      nmat=np+1
      if(nmat.gt.npm) then
       write(6,*) 'ERROR, recompile with npm>= ',nmat
       stop
      endif

      write(6,*) 'nmat read=',nmat

      rewind(12) 

      ngen=10000000

      wbin=0
      do kk=1,nmat
       ebin2(kk)=0.d0
       ebin(kk)=0.d0
      enddo
 
      ibin=0

      do kt=1,ngen
       read(12,end=100) (ebuf(kk),kk=1,nmat)
       ibin=ibin+1
       if(ibin.ge.ibinit) then
        wbin=wbin+1.d0
        do kk=1,np
         ebin(kk)=ebin(kk)+ebuf(kk)
         ebin2(kk)=ebin2(kk)+ebuf(kk)**2
        enddo
        write(21,123) ibin,(ebuf(kk),kk=1,np)
       endif 
      enddo 

100   continue

      nmis=ibin-ibinit+1
      write(6,*) ' number of measures done =',nmis

      do kk=1,np
       ebin(kk)=ebin(kk)/wbin
       ebin2(kk)=dsqrt(dabs(ebin2(kk)/wbin-ebin(kk)**2))
       ebin2(kk)=ebin2(kk)/dsqrt(dfloat(nmis))
      enddo
      write(20,124) (ebin(kk),ebin2(kk),kk=1,np)
      write(20,*) ' with no error bars '
      write(20,124) (ebin(kk),kk=1,np)

123   format(i7,1000e16.7)        
124   format(1000e16.7)        

      close(13)
      close(20) 
      close(21) 

      stop
      end
