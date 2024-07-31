        program bootback 
        implicit real*8(a-h,o-z)
        parameter(nh=500,nbinm=10000000,nm=1)
        dimension  w(nbinm),e(nbinm),es(nbinm)
         iseed=12375731
         call rand_init(iseed)

c       y =  a*x*2+ b*x + c
c       legge z,y da for011.dat 
        eps=1d-6
        open(unit=31,form='formatted',status='old',file='fort.31')
        open(unit=32,form='formatted',status='old',file='fort.32')

c        rewind(20) 
	n=1
C	lettura dati


        nbin=0
        dowhile(nbin.lt.nbinm) 
        nbin=nbin+1
5       read(31,*,end=500) e(nbin),w(nbin)
        read(32,*,end=500) es(nbin),wpip
        if(wpip.ne.w(nbin)) then
        write(6,*) ' warning  the measures are not correlated !!!! '
        endif
        enddo
500        continue


       if(nbin.ne.nbinm) then
       nbin=nbin-1
       else
       write(6,*) ' Warning maximum number of bins exceeded !!!'
       endif




       write(6,*) ' number of bins read =',nbin 


        nmis=1000

c       initialization opara e opars
        err=0.d0
        eta=0.d0

        eav=0.d0
        esav=0.d0 

         
       do kmain=1,nmis

c
        wt=0.d0
        et=0.d0
        ets=0.d0

c       bootstrap estimate sample 
        do k=1,nbin
        j=drand1()*nbin+1 
c       j=k 
         wt=wt+w(j)
         et=et+w(j)*e(j)
         ets=ets+w(j)*es(j)

        enddo

         et=et/wt
         ets=ets/wt       


c       now average correlation function on the given bootstrap


        givenmeas=ets-et**2

c       write(6,*) ' givenmeas =',givenmeas 

c       givenmeas=opart(i,ii)

        eta=eta+givenmeas
        err=err+givenmeas**2

        eav=eav+et
        esav=esav+et**2

       enddo

      eta=eta/nmis
      err=err/nmis

      esav=esav/nmis
      eav=eav/nmis

      err=dsqrt(err-eta**2)
      esav=dsqrt(esav-eav**2)


        write(6,*)  ' Energy =',eav,esav
        write(6,*)  ' Variance square =',eta,err

        stop
	end

