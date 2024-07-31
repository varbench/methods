            program erread
            implicit none
            integer :: maxk,k,kt,i,j,ng,nbin &
      ,iskip,maxj,ibinit,icount,ibin,nmis,lbin,iskipr,cont
            integer,parameter:: npm=1000,ngen=20000000
       real*8,dimension(:),allocatable:: ebuf,wbuf,ek,wk &
       ,ebin,wbin,ebin2
            real*8,dimension(:),allocatable:: wbufm,wtot
            real*8 :: eskip(npm),wsk
            
          write(6,*) ' max k corrections, bin lenght, ibinit,iskip ?' 
           read(5,*) maxk,lbin,ibinit,iskipr
       
      allocate(ebuf(ngen),wbuf(ngen),ek(0:ngen),wk(0:ngen),ebin(0:ngen) &
     ,wbin(0:ngen),ebin2(0:ngen),wbufm(ngen),wtot(ngen)) 
      if(iskipr.ge.0) then 
      open(unit=12,file='fort.12',form='unformatted',status='old')
      iskip=iskipr
      else
      iskip=-iskipr-1
      open(unit=12,file='fort.12.fn',form='unformatted',status='old')
      endif


      open(unit=20,file='fort.30',form='formatted',status='unknown')
      open(unit=21,file='fort.31',form='formatted',status='unknown')
            rewind(12) 


           if(iskip.gt.npm-1) then
           write(6,*) '!!! ERRORE recompile with npm> ',iskip
           stop
           endif

!          iskip=0   --> energy
!          iskip=1   --> variance 
!          iskip=2   --> order parameter   
           if(iskip.eq.0) then 
           write(20,*) ' Calculation gs energy '
           endif

           do i=0,ngen
           ek(i)=0.d0
           wk(i)=0.d0
           ebin(i)=0.d0
           wbin(i)=0.d0
           ebin2(i)=0.d0
           enddo

          
           icount=0
           ibin=0
           
           kt=1         
              do i=1,ngen
              wbufm(i)=wbuf(i)
              enddo

              maxj=ngen
              cont=0

              do j=1,maxj
              if(iskip.ne.0) then
       read(12,end=100) wbuf(j),wtot(j),(eskip(k),k=1,iskip),ebuf(j)
              else
              read(12,end=100) wbuf(j),wtot(j),ebuf(j)
              endif
              cont=cont+1
              enddo
100        continue
    
             maxj=cont  
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
        
              enddo   ! end j

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
                 wsk=wsk*wbufm(j-k+ngen)                     
                 endif
                endif 
                ek(k)=ek(k)+ebuf(j)*wsk
                wk(k)=wk(k)+wsk
                enddo      ! end do k 

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
             enddo 
             endif 
             
!            calculation error bars

             nmis=ibin-ibinit+1
             write(6,*) ' number of measures done =',nmis
             do i=0,maxk
             ebin(i)=ebin(i)/wbin(i)
             ebin2(i)=dsqrt(dabs(ebin2(i)/wbin(i)-ebin(i)**2))
             ebin2(i)=ebin2(i)/dsqrt(dfloat(nmis))
             enddo
          write(20,*) ' Independent bins ',nmis,'of lenght ',lbin
          write(20,*) 
           write(20,*) ' Energy , error, # of bias  correcting factor '
           do i=0,maxk
           write(20,*) ebin(i),ebin2(i),i
           enddo
           close(11) 
           close(12)
           close(20) 
        
              stop
              end
      
