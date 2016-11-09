c=======================================================================

	program main

	implicit none

	integer n,param1
	character*80  file,label1
        character*80  arg
        integer plot
        integer nargs

c-----------------------------------------------------------------------
c
c       Check the command line argument 
        plot=1 
        nargs = iargc()
        if( nargs.gt.0 ) then
                call getarg( 1, arg )
                if( arg.eq.'-n' ) then
                        plot=0
                endif
        endif


	write(*,*) 
	open(unit=8,file='histogram.inp',status='unknown')
	call options(file,param1,label1,n)
	call driver(file,param1,label1,n,plot)
	close(8)
	write(*,*) 
	write(*,*) 'Input commands saved in histogram.inp'
	write(*,*) 

	end

c=======================================================================

	subroutine options(file,param1,label1,n)
	
	implicit none
	
	integer n,param1
	character*80  file,label1

c-----------------------------------------------------------------------
	
	write(*,'(a,$)') ' Input filename (def=bayes.dat): '
	read(*,'(a)') file
	if( file.eq.'') then 
	  file = 'bayes.dat'
	endif
	write(8,'(a80)') file

	write(*,'(a,$)') ' Input Column : '
	read(*,*) param1
	write(8,*) param1
	write(*,'(a,$)') ' Input x axis label: '
	read(*,'(a)') label1
	write(8,'(a80)') label1

	call nlines(file,n)
	

	return
	end
	
c=======================================================================
c Count the number of lines in the bayes.dat file
c=======================================================================
        subroutine nlines(filename,lines)
c       
        integer lines
        character*80 filename
        character*80 record
c
        open(unit=9,file=filename,form='formatted',status='old')
c       
        lines = 0
        do i=1,1000000
          read(9,'(a)',end=2) record 
          lines = lines + 1
        enddo
c
2       close(9)
        return
        end

c=======================================================================

	subroutine driver(file,param1,label1,n,plot)
	
	implicit none
	
        include 'histogram.inc'
        
	integer  n,n2,param1
        integer  plot
	character*80  file,ans,label1,label2,title
	character record* (REC_SIZE)

	real x(n),fld(NFS)
	real xbar,xerr
	real percent
	real lower67,upper67
	real lower95,upper95
	real lower99,upper99
	integer nbin,large,nx2,current
	parameter(large=1000,nx2=256)
	real hist(large),cell,tr(2),shist(256)
	real xmin,xmax,ymin,ymax,x1,x2,y1,y2,xold,yold
	real frac(3),levels(3),sell,xp(4),yp(4),fwhm
	real u(nx2-1),f(nx2-1),a,b,norm,peak
	integer i,j,k,ipeak,idist
        

c-----------------------------------------------------------------------

	label2 = 'Probability density'
	title = ' '
        if( plot.eq.1) then
	  call pgbeg(0,'/xs',1,1)
        else
	  call pgbeg(0,'/null',1,1)
        endif
	call pgask(.false.)
	call pgpage
	call pgpap(0.0,0.707106781)
	call pgsch(1.5)
	call pgslw(3)

c Read in data:

	open(unit=9,file=file,form='formatted',status='old')
       
	n2 = 1
	do i=1,n
	  read(9,'(a)') record
	  if( record(1:1).ne.'#' ) then
            read(record,*,ERR=8,END=8) fld
8           CONTINUE
            x(n2) = fld(param1)
	    n2 = n2 + 1
	  endif
	enddo
	close(9)
	n = n2 - 1

c Do stats:

	call sort(x,n)
	
	xmin = x(n)
	xmax = x(1)
	i = 0.165*n

	call stats(n,x,xbar,xerr)

	write(*,*) 'Statistics:'
	write(*,*) '  Ndata = ',n
	write(*,*) '  Mean = ',xbar
	write(*,*) '  Std. dev. = ',xerr
	write(*,*) '  Median = ',x(n/2)
	write(*,*) '  16.5 percentile = ',x(n-i)
	write(*,*) '  83.5 percentile = ',x(i)
	write(*,*)

c Bin:

	write(*,*) 'x limits are (',xmin,' - ',xmax,')'
        write(*,'(a,$)') ' Change? (y or n, def=n): '
	read(*,'(a1)') ans
	write(8,'(a1)') ans
	if (ans.eq.'y') then
	  write(*,'(a,$)') '   Enter new limits (xmin,xmax): '
	  read(*,*) xmin,xmax
	  write(8,*) xmin,xmax 
	endif
	write(*,*)
	 
10	write(*,'(a,$)') ' Enter no. of bins: '
	read(*,*) nbin
c	if (mod(nbin,2).ne.0) then
c        write(*,*) '  Nbin must be even...'
c	  goto 10
c	endif  
	
        write(*,*) 'Making histogram...'
	
	cell = (xmax-xmin)/float(nbin-1)
	tr(1) = xmin - cell
	tr(2) = cell
	
	write(*,*) 
	do i=1,nbin
	  hist(i) = 0.0
	enddo
	
	call bin(n,x,nbin,hist,tr)
	
	ymin = 0.0
	ymax = -1e32
	do i=1,nbin
	  ymax = max(ymax,hist(i))
	enddo
	ymax = ymax*1.1
	
c check sensibleness of contours

	call findlevels(nbin,cell,levels,hist)
	do k=1,3
	  frac(k) = 0.0
	enddo
        do i=1,nbin
	  do k=1,3
	    if (hist(i).ge.levels(k)) then
	      frac(k) = frac(k) + hist(i)*cell
	    endif
	  enddo
	enddo
	do k=1,3
	  frac(k) = 100.0*frac(k)
	enddo
 	write(*,'(a,f5.1,a)') 
     &'    68% contour contains ',frac(1),'% of samples,'
 	write(*,'(a,f5.1,a)') 
     &'    95% contour contains ',frac(2),'% of samples,'
 	write(*,'(a,f5.1,a)') 
     &'    99% contour contains ',frac(3),'% of samples,'
 	write(*,*) 
	frac(1) = frac(1)-68 + frac(2)-95 + frac(3)-99
	write(*,'(a,f5.1)')
     &'    error : ',frac(1)

	call pgpage
	call pgvport(0.15,0.85,0.15,0.85)
	call pgwindow(xmin,xmax,ymin,ymax)
	xold = xmin
	yold = ymax
	do i=1,nbin
	  x1 = tr(1) + i*tr(2)
	  x2 = x1 + cell
	  y1 = 0.0
	  y2 = hist(i)
	  call pgsci(1)
	  if (hist(i).ge.levels(3)) call pgsci(7)
	  if (hist(i).ge.levels(2)) call pgsci(8)
	  if (hist(i).ge.levels(1)) call pgsci(2)
	  call pgrect(x1,x2,y1,y2)
	  call pgsci(1)
	  call pgmove(xold,yold)
	  call pgdraw(x1,y2)
	  call pgdraw(x2,y2)
	  xold = x2
	  yold = y2
	enddo
	call pgsci(1)
	call pgbox('BCTSN',0.0,0,'BCTSN',0.0,0)
	call pglab(label1,label2,title)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c Offer to rebin the data:
	
12	write(*,'(a,$)') ' Re-bin? (y or n): '
	read(*,'(a1)') ans
	if (ans.eq.'y') then 
	  goto 10
	endif
	write(8,*) nbin
	write(8,'(a1)') 'n'
	write(*,*)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c Offer to overplot some other distribution:

14	idist = 0
	write(*,'(a,$)') ' Overlay standard distribution? (y or n): '
	read(*,'(a1)') ans
	if (ans.ne.'y') then
	  goto 18
	endif

	do i=1,nx2-1
	  u(i) = xmin + (xmax-xmin)*float(i-1)/float(nx2-2)
	enddo
	
	write(*,'(a)') 
16	write(*,'(a)') ' Overlay options: '
	write(*,'(a)') '   1) uniform'
	write(*,'(a)') '   2) Jeffreys'
	write(*,'(a)') '   3) Gaussian'
	write(*,'(a)') '   4) Lognormal'
	write(*,'(a)') '   5) Exponential'
	write(*,'(a)') '   6) Cauchy/Lorentzian'
	write(*,'(a,$)') ' Enter option number (1..6): '
	read(*,*) idist

	if (idist.eq.1) then
	  write(*,'(a,$)') ' Enter range limits: '
	  read(*,*) a,b
	  do i=1,nx2-1
	    if (u(i).ge.a.and.u(i).lt.b) then
	      f(i) = 1.0/(b-a)
	    else
	      f(i) = 0.0
	    endif
	  enddo
	elseif (idist.eq.2) then
	  write(*,'(a,$)') ' Enter range limits: '
	  read(*,*) a,b
	  norm = 1.0/log(b/a)
	  do i=1,nx2-1
	    if (u(i).ge.a.and.u(i).lt.b) then
	      f(i) = norm/u(i)
	    else
	      f(i) = 0.0
	    endif
	  enddo
	elseif (idist.eq.3) then
	  write(*,'(a,$)') ' Enter mean and sigma: '
	  read(*,*) a,b
	  norm = 1.0/sqrt(6.283185307*b*b)
	  do i=1,nx2-1
	    f(i) = norm*exp(-(u(i)-a)*(u(i)-a)/(2.0*b*b))
	  enddo
	elseif (idist.eq.4) then
	  write(*,'(a,$)') ' Enter peak posn and width: '
	  read(*,*) a,b
	  norm = exp(-0.5*b*b)/(a*sqrt(6.283185307*b*b))
	  do i=1,nx2-1
	    f(i) = norm*exp(-((log(u(i))-log(a))**2.0)/(2.0*b*b))
	  enddo
	elseif (idist.eq.5) then
	  write(*,'(a,$)') ' Enter mean: '
	  read(*,*) a
	  do i=1,nx2-1
	    f(i) = exp(-u(i)/a)/a
	  enddo
	elseif (idist.eq.6) then
	  write(*,'(a,$)') ' Enter mean and half width: '
	  read(*,*) a,b
	  norm = b/3.141592654
	  do i=1,nx2-1
	    f(i) = norm/(b*b - (u(i)-a)*(u(i)-a))
	  enddo
	else
	  goto 16
	endif
      
	call pgeras
	call pgvport(0.15,0.85,0.15,0.85)
	call pgwindow(xmin,xmax,ymin,ymax)
	xold = xmin
	yold = ymax
	do i=1,nbin
	  x1 = tr(1) + i*tr(2)
	  x2 = x1 + cell
	  y1 = 0.0
	  y2 = hist(i)
	  call pgsci(1)
	  if (hist(i).ge.levels(3)) call pgsci(7)
	  if (hist(i).ge.levels(2)) call pgsci(8)
	  if (hist(i).ge.levels(1)) call pgsci(2)
	  call pgrect(x1,x2,y1,y2)
	  call pgsci(1)
	  call pgmove(xold,yold)
	  call pgdraw(x1,y2)
	  call pgdraw(x2,y2)
	  xold = x2
	  yold = y2
	enddo
	call pgsci(5)
	call pgline(nx2-1,u,f)
	call pgsci(1)
	call pgbox('BCTSN',0.0,0,'BCTSN',0.0,0)
	call pglab(label1,label2,title)
            
17	write(*,'(a,$)') ' Try again? (y or n): '
	read(*,'(a1)') ans
	if (ans.eq.'n') then
	  write(8,'(a1)') 'y'
	  write(8,*) idist
	  if (idist.eq.5) then
	    write(8,*) a
	  else
	    write(8,*) a,b
	  endif
	  write(8,'(a1)') 'n'
	  call pgend
	elseif (ans.eq.'y') then
	  goto 16
	else 
	  goto 17
	endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c Offer to replot as postscript:
	
18	write(*,*)
	write(*,'(a,$)') ' Plot postcript? (y or n): '
	read(*,'(a1)') ans
	write(8,'(a1)') ans 
	if (ans.eq.'y') then 
	  call pgbeg(0,'/vcps',1,1)
	  call pgask(.false.)
	  call pgpap(0.0,0.707106781)
	  call pgpage
	  call pgsch(1.5)
	  call pgslw(3)
	  call pgvport(0.15,0.85,0.15,0.85)
	  call pgwindow(xmin,xmax,ymin,ymax)
	  xold = xmin
	  yold = ymax
	  do i=1,nbin
	    x1 = tr(1) + i*tr(2)
	    x2 = x1 + cell
	    y1 = 0.0
	    y2 = hist(i)
	    call pgsci(7)	      
	    if (hist(i).ge.levels(3)) call pgsci(7)	     
	    if (hist(i).ge.levels(2)) call pgsci(8)	      
	    if (hist(i).ge.levels(1)) call pgsci(2)	      
	    call pgrect(x1,x2,y1,y2)
	    call pgsci(1)
	    call pgmove(xold,yold)
	    call pgdraw(x1,y2)
	    call pgdraw(x2,y2)
	    xold = x2
	    yold = y2
	  enddo
	  if (idist.gt.0) then
	    call pgsci(5)
	    call pgline(nx2-1,u,f)
	  endif
	endif
	call pgsci(1)
	write(*,'(a,$)') ' Plot points? (y or n): '
	read(*,'(a1)') ans
	write(8,'(a1)') ans 
	if (ans.eq.'y') then 	
	  do i=1,n
	    call pgmove(x(i),ymin)
	    call pgdraw(x(i),ymin+(ymax-ymin)/25.0)
	  enddo  
	endif
	write(*,*) 
	write(*,'(a,$)') ' Plot reference point? (y or n): '
	read(*,'(a1)') ans
	write(8,'(a1)') ans 
	if (ans.eq.'y') then 	
	  write(*,'(a,$)') '    Enter value: '
	  read(*,*) x1
	  write(8,*) x1 
	  call pgsls(4)
	  call pgslw(3)
	  call pgmove(x1,ymin)
	  call pgdraw(x1,ymax)
	  call pgsls(1)
	  call pgslw(1)
	endif
	call pgbox('BCTSN',0.0,0,'BCTSN',0.0,0)
	call pglab(label1,label2,title)
	call pgend

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c Write out overlay data:

	if (idist.gt.0) then
        open(unit=9,file='overlay.txt',status='unknown')
        do i=1,nx2-1
          write(9,*) u(i),f(i)
        enddo
        close(9)
	endif

c Write out histogram data:

        open(unit=9,file='histogram.txt',status='unknown')
        do i=1,nbin-1
          u(i) = tr(1) + i*tr(2)
          f(i) = hist(i)
          write(9,*) u(i),f(i)
        enddo
        close(9)

c-----------------------------------------------------------------------

	write(*,*) 
	write(*,'(a,$)') ' Smooth? (y or n): '
	read(*,'(a1)') ans
	write(8,'(a1)') ans 
	if (ans.eq.'y') then 
	  call pgbeg(0,'/xs',1,1)
          call pgask(.false.)
          call pgpap(0.0,0.707106781)
	  call pgsch(1.5)
	  
20	  write(*,'(a,$)') '   Enter smoothing FWHM: '
	  read(*,*) fwhm
	  
	  call smoothhist(nbin,large,hist,xmin,xmax,cell,fwhm,shist)
	  sell = cell*nbin/float(nx2)
	  call findlevels(nx2,sell,levels,shist)
	  call findpeak(nx2,shist,ipeak)
	  
	  ymin = 0.0
	  ymax = -1e32
	  do i=1,nx2
	    ymax = max(ymax,shist(i))
	  enddo
	  ymax = ymax*1.1
	  tr(1) = xmin - sell
	  tr(2) = sell
	  peak = tr(1)+ipeak*tr(2)+0.5*sell
   	  write(*,*) '   Peak value = ',peak
          current = 1
	  do i=1,nx2
	    if     (shist(i).ge.levels(3).and.current.eq.1) then
	      lower99 = tr(1)+i*tr(2)+0.5*sell
	      current = 2
	    elseif (shist(i).ge.levels(2).and.current.eq.2) then
	      lower95 = tr(1)+i*tr(2)+0.5*sell
	      current = 3
	    elseif (shist(i).ge.levels(1).and.current.eq.3) then
	      lower67 = tr(1)+i*tr(2)+0.5*sell
	      current = 4
	    elseif (shist(i).le.levels(1).and.current.eq.4) then
	      upper67 = tr(1)+i*tr(2)+0.5*sell
	      current = 5
	    elseif (shist(i).le.levels(2).and.current.eq.5) then
	      upper95 = tr(1)+i*tr(2)+0.5*sell
	      current = 6
	    elseif (shist(i).le.levels(3).and.current.eq.6) then
	      upper99 = tr(1)+i*tr(2)+0.5*sell
	      current = 7
	    endif
	    if (current.eq.7) goto 30
	  enddo
30	  write(*,*) '  68% interval = ',lower67,' ... ',upper67
	  write(*,*) '  95% interval = ',lower95,' ... ',upper95
	  write(*,*) '  99% interval = ',lower99,' ... ',upper99
	  write(*,*) ' '
	  write(*,*) '  68% error bars: ',lower67-peak,', +',upper67-peak

          call pgpage
	  call pgvport(0.15,0.85,0.15,0.85)
          call pgwindow(xmin,xmax,ymin,ymax)
	  do i=1,nx2-1
	    xp(1) = tr(1) + i*tr(2)
	    xp(2) = xp(1) + sell
	    xp(3) = xp(2)
	    xp(4) = xp(1)
	    yp(1) = 0.0
	    yp(2) = 0.0
	    yp(3) = shist(i+1)
	    yp(4) = shist(i)
	    call pgsci(1)
	    if (shist(i).ge.levels(3)) call pgsci(7)
	    if (shist(i).ge.levels(2)) call pgsci(8)
	    if (shist(i).ge.levels(1)) call pgsci(2)
	    call pgpoly(4,xp,yp)
	    call pgsci(1)
	    call pgmove(xp(3),yp(3))
	    call pgdraw(xp(4),yp(4))
          enddo
	  call pgsci(1)
	  xp(1) = tr(1) + ipeak*tr(2)
	  call pgarro(xp(1)+0.5*sell,0.1*ymax,xp(1)+0.5*sell,0.0)
	  call pgbox('BCTSN',0.0,0,'BCTSN',0.0,0)
          call pglab(label1,label2,title)
	  
	  write(*,'(a,$)') ' Re-smooth? (y or n): '
	  read(*,'(a1)') ans
	  if (ans.eq.'y') then
	    goto 20
	  endif
	  write(8,*) fwhm 
	  write(8,'(a1)') 'n' 
	  call pgend

c Write out smooth curve data:

          open(unit=9,file='smooth_histogram.txt',status='unknown')
c 	  write(9,*) nx2-1
c 	  write(9,*) 
	  do i=1,nx2-1
	    u(i) = tr(1) + i*tr(2)
	    f(i) = shist(i)
	    write(9,*) u(i),f(i)
          enddo
	  close(9)
	  
	  write(*,'(a,$)') ' Plot postcript? (y or n): '
	  read(*,'(a1)') ans
	  write(8,'(a1)') ans
	  if (ans.eq.'y') then 	
	  call pgbeg(0,'smooth_histogram.ps/vcps',1,1)
          call pgask(.false.)
          call pgpap(0.0,0.707106781)
          call pgpage
	  call pgsch(1.5)
	  call pgvport(0.15,0.85,0.15,0.85)
          call pgwindow(xmin,xmax,ymin,ymax)
	  do i=1,nx2-1
	    xp(1) = tr(1) + i*tr(2)
	    xp(2) = xp(1) + sell
	    xp(3) = xp(2)
	    xp(4) = xp(1)
	    yp(1) = 0.0
	    yp(2) = 0.0
	    yp(3) = shist(i+1)
	    yp(4) = shist(i)
	    call pgsci(7)
	    if (shist(i).ge.levels(3)) call pgsci(10)
	    if (shist(i).ge.levels(2)) call pgsci(11)
	    if (shist(i).ge.levels(1)) call pgsci(12)
	    call pgpoly(4,xp,yp)
	    call pgsci(1)
	    call pgmove(xp(3),yp(3))
	    call pgdraw(xp(4),yp(4))
          enddo
	  call pgsci(1)
	  xp(1) = tr(1) + ipeak*tr(2)
	  call pgarro(xp(1)+0.5*sell,0.1*ymax,xp(1)+0.5*sell,0.0)
	  write(*,'(a,$)') ' Plot points? (y or n): '
	  read(*,'(a1)') ans
	  write(8,'(a1)') ans
	  if (ans.eq.'y') then 	
	    do i=1,n
	      call pgmove(x(i),ymin)
	      call pgdraw(x(i),ymin+(ymax-ymin)/25.0)
	    enddo  
	  endif
	  write(*,*) 
	  write(*,'(a,$)') ' Plot reference point? (y or n): '
	  read(*,'(a1)') ans
	  write(8,'(a1)') ans
	  if (ans.eq.'y') then
	    write(*,'(a,$)') '    Enter value: '
	    read(*,*) x1
	    write(8,*) x1
	    call pgsls(4)
	    call pgslw(3)
	    call pgmove(x1,ymin)
	    call pgdraw(x1,ymax)
	    call pgsls(1)
	    call pgslw(1)
	  endif
	  call pgbox('BCTSN',0.0,0,'BCTSN',0.0,0)
          call pglab(label1,label2,title)
	  call pgend
	  endif
	  
	endif

	return
	end
	
c=======================================================================

        subroutine stats(n,x,xbar,xerr)

        implicit none
       
        integer n
        real x(n),xbar,xerr

        integer i
        real*8 sum,sumsq,nn

c-----------------------------------------------------------------------

        nn = float(n)

        sum = 0.0
        sumsq = 0.0     
        do i=1,n
          sum = sum + x(i)
c          sumsq = sumsq + x(i)*x(i)
        enddo
        xbar = sum/nn  

	do i=1,n
          sumsq = sumsq + (x(i)-xbar)*(x(i)-xbar)
        enddo

	xerr = sqrt(sumsq/(nn-1.0))

c        xerr = sumsq - sum*sum/nn
c        xerr = sqrt(xerr/(nn-1.0))

        return
        end
       
c=======================================================================
c bin data, and return hist array as an estimate of probability density

	subroutine bin(n,x,nbin,hist,tr)
	
	implicit none
	
	integer n,nbin
	real x(n),hist(1000),tr(2)

	real x1,x2
	integer i,j

	do i=1,nbin         
	  x1 = tr(1) + i*tr(2)
	  x2 = x1 + tr(2)
	  hist(i) = 0.0
	  do j=1,n
	    if (x(j).gt.x1.and.x(j).le.x2) then
	      hist(i) = hist(i) + 1.0
	    endif
	  enddo
	  hist(i) = hist(i)/(float(n)*tr(2))
	enddo

      return
      end
       
c=======================================================================

	include 'sort.f'

c=======================================================================

	subroutine findlevels(n,cell,levels,z)
	
	implicit none
	
	integer n
	real z(*),copy(n),levels(3),cell
	
	real zmax,sum,sum2
	integer i,count,current
	
	sum = 0.0
	do i=1,n	
	  copy(i) = z(i)
	  sum = sum + copy(i)*cell
	enddo
	do i=1,n	
	  copy(i) = copy(i)/sum
	enddo
	
	call sort(copy,n)
	
      current = 1
	sum2 = 0.0
	do i=1,n
	  sum2 = sum2 + copy(i)*cell
	  if (sum2.ge.0.68.and.current.eq.1) then
	    levels(1) = copy(i)*sum
	    current = 2
	  elseif (sum2.ge.0.95.and.current.eq.2) then
	    levels(2) = copy(i)*sum
	    current = 3
	  elseif (sum2.ge.0.99.and.current.eq.3) then
	    levels(3) = copy(i)*sum
	    current = 4
	  endif
	  if (current.eq.4) goto 999
	enddo
	
999	return
	end
	
c=======================================================================
      
	subroutine smoothhist(nx,large,array,xmin,xmax,cell,fwhm,
     &                      smoothcrop)
	
	implicit none
	
	integer nx,large
	real array(large),array2(2*nx),xmin,xmax,cell,fwhm
	integer nx2,nx3,iw,npix,i
	parameter(nx2=256,nx3=2*nx2)
	real smootharray(nx3),smoothcrop(nx2),sum,sell
	
	npix = int((fwhm*nx2)/(cell*nx)) + 1
	sell = cell*nx/float(nx2)
	
c Fill working array, pad with zeroes

	do i=1,nx/2
	  array2(i) = 0.0
	enddo
	do i=nx/2+1,3*nx/2
	  array2(i) = array(i-nx/2)
	enddo
	do i=3*nx/2+1,2*nx
	  array2(i) = 0.0
	enddo
	
	iw = 0
	call smooth1D(array2,2*nx,smootharray,nx3,npix,iw)
	
	sum = 0.0
	do i=1,nx3
 	  sum = sum + smootharray(i)*sell
	enddo
	do i=1,nx3
c 	  smootharray(i) = smootharray(i)*sell/sum
	  smootharray(i) = smootharray(i)/sum
	enddo
	
c Crop smoothed array

	do i=1,nx2
	  smoothcrop(i) = smootharray(i+nx2/2)
	enddo
	
	return
	end
	
c=======================================================================

	subroutine findpeak(n,z,ipeak)
	
	integer n,ipeak
	real z(*),copy(n)
	
	real zmax,sum
	integer i,count,current
	
	zmax = 0.0
	do i=1,n	
	  if (z(i).gt.zmax) then
	   zmax = z(i)
	   ipeak = i
	  endif
	enddo

	return
	end

c=======================================================================

	include 'smooth1D.f'

c=======================================================================

