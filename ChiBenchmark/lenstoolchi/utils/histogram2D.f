c=======================================================================
c
	program main
c
	implicit none
c	
	integer npts,n,param1,param2
	character*80  file1,label1,label2,title
        character*80  arg
        integer plot
        integer nargs
c
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

	open(unit=8,file='histogram2d.inp',status='unknown')
        call options(file1,param1,param2,npts,n,label1,label2,title)
	call driver(file1,param1,param2,npts,n,label1,label2,title,plot)
	close(8)
	write(*,*) 
	write(*,*) 'Input commands saved in histogram2d.inp'
	write(*,*) 
	
c
	end
c		
c=======================================================================
c
	subroutine options(file1,param1,param2,n,nbin,label1,label2,title)
c	
	implicit none
c	
	integer n,nbin,param1,param2
	character*80 file1,label1,label2,title
c
c-----------------------------------------------------------------------
c	
	write(*,'(a,$)') ' Input filename (def=bayes.fits) : '
	read(*,'(a)') file1
	if( file1.eq.'' ) then
	  file1 = 'bayes.fits'
	endif
	write(8,'(a80)') file1
	write(*,'(a,$)') ' Input Column 1 : '
	read(*,*) param1
	write(8,*) param1
	write(*,'(a,$)') '   Input plot label 1: '
	read(*,'(a)') label1
	write(8,'(a80)') label1
	write(*,'(a,$)') ' Input Column 2 : '
	read(*,*) param2
	write(8,*) param2
	write(*,'(a,$)') '   Input plot label 2: '
	read(*,'(a)') label2
	write(8,'(a80)') label2
	write(*,'(a,$)') ' Input title for plot: '
	read(*,'(a)') title
	write(8,'(a80)') title

        call nlines(file1,n)

	write(*,*) 'No. of lines = ',n,'...'
	write(*,'(a,$)') ' ...input no. of pixels: '
	read(*,*) nbin
	write(8,*) nbin

	write(*,*)
	
	return
	end
	
c=======================================================================
c
        subroutine nlines(filename,lines)
c       
        integer lines
        integer status,unit,readwrite,blocksize,hdutype
        character*80 filename

        status = 0
        call ftgiou(unit, status)
        readwrite = 0
        call ftopen(unit,filename,readwrite,blocksize,status)
        call ftmahd(unit,2,hdutype,status)
        call ftgnrw(unit, lines, status)
        call ftclos(unit, status)
        call ftfiou(unit, status)
        if (status .gt. 0)call printerror(status)
        return
        end
c
c *************************************************************************
        subroutine printerror(status)

c  This subroutine prints out the descriptive text corresponding to the
c  error status value and prints out the contents of the internal
c  error message stack generated by FITSIO whenever an error occurs.

        integer status
        character errtext*30,errmessage*80

c  Check if status is OK (no error); if so, simply return
        if (status .le. 0)return

c  The FTGERR subroutine returns a descriptive 30-character text string that
c  corresponds to the integer error status number.  A complete list of all
c  the error numbers can be found in the back of the FITSIO User s Guide.
        call ftgerr(status,errtext)
        print *,'FITSIO Error Status =',status,': ',errtext

c  FITSIO usually generates an internal stack of error messages whenever
c  an error occurs.  These messages provide much more information on the
c  cause of the problem than can be provided by the single integer error
c  status value.  The FTGMSG subroutine retrieves the oldest message from
c  the stack and shifts any remaining messages on the stack down one
c  position.  FTGMSG is called repeatedly until a blank message is
c  returned, which indicates that the stack is empty.  Each error message
c  may be up to 80 characters in length.  Another subroutine, called
c  FTCMSG, is available to simply clear the whole error message stack in
c  cases where one is not interested in the contents.
        call ftgmsg(errmessage)
        do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
        end do
        end

c=======================================================================
	subroutine driver(file1,param1,param2,npts,n,label1,label2,title,plot)
	
	implicit none

        include 'histogram.inc'
	
	integer  n,npts,param1,param2
        integer  plot
	character record* (REC_SIZE)
	character*80  file1
	character*80 label1,label2,title,ans

	real fwhm
	real x(npts),y(npts),fld(NFS)
	real pt1,pt2,err1,err2
	integer n2
	real hist(n,n),hist2(2*n,2*n),shist(n,n)
	real xmin,xmax,ymin,ymax,xcentre,ycentre,xcell,ycell,tr(6)
	real x0,y0,u,zmin,zmax,lmin,lmax,sum,frac(3),levels(3)
	integer i,j,k,nc,iw
	parameter (nc=5)
	real c(nc),r(nc),g(nc),b(nc)
	real smin,smax
	real minerror, bestfwhm
        integer status,unit,readwrite,blocksize,hdutype
        integer felem,frow
        real nulle
        logical anynull

c-----------------------------------------------------------------------


c Read in data:
	write(*,*) 'Reading in data...'
	write(*,*) 

        status = 0
        call ftgiou(unit,status)
        readwrite = 0
        call ftopen(unit,file1,readwrite,blocksize,status)
        call ftmahd(unit,2,hdutype,status)
        frow = 1
        felem = 1
        nulle = 0.
        call ftgcve(unit,param1,frow,felem,npts,nulle,x,anynull,status)
        call ftgcve(unit,param2,frow,felem,npts,nulle,y,anynull,status)
        call ftclos(unit, status)
        call ftfiou(unit, status)

	xmin = 1e32
	xmax = -1e32
	ymin = 1e32
	ymax = -1e32
        do i=1,npts
	    xmin = min(x(i),xmin)
	    xmax = max(x(i),xmax)
	    ymin = min(y(i),ymin)
	    ymax = max(y(i),ymax)
        enddo
        write(*,*) npts,' lines read'

	write(*,*) 'x limits are (',xmin,' - ',xmax,')' 
	write(*,'(a,$)') ' Change? (y or n, def=n): '
	read(*,'(a1)') ans
	write(8,'(a1)') ans
	if (ans.eq.'y') then 
	  write(*,'(a,$)') '   Enter new limits (xmin,xmax): '
	  read(*,*) xmin,xmax
	  write(8,*) xmin,xmax 
	  write(*,*) 
	endif
	write(*,*) 'y limits are (',ymin,' - ',ymax,')' 
	write(*,'(a,$)') ' Change? (y or n, def=n): '
	read(*,'(a1)') ans
	write(8,'(a1)') ans
	if (ans.eq.'y') then
	  write(*,'(a,$)') '   Enter new limits (ymin,ymax): '
	  read(*,*) ymin,ymax
	  write(8,*) ymin,ymax 
	endif
	write(*,*) 

	xcell = (xmax-xmin)/float(n-1)
	ycell = (ymax-ymin)/float(n-1)
	tr(1) = xmin - xcell
	tr(2) = xcell
	tr(3) = 0.0
	tr(4) = ymin - ycell
	tr(5) = 0.0
	tr(6) = ycell

	n2 = 2*n
	write(*,*) 'Making histogram...'
	write(*,*) 
	do i=1,n
        do j=1,n
	    hist(i,j) = 0.0
	  enddo
	enddo
	do i=1,n2
        do j=1,n2
	    hist2(i,j) = 0.0
	  enddo
	enddo
	
        call bin(n,n,npts,x,y,hist,xcell,ycell,tr)
	call pad(n,n,hist,n2,n2,hist2)
	
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - -

 10	write(*,'(a,$)') ' Look for the best FWHM? (y or n) : '
	read(*,'(a1)') ans
	if( ans.eq.'y' ) then
	  goto 11
	else
	  write(*,'(a,$)') ' Enter FWHM (pixels): '
	  read(*,*) fwhm
	  goto 15
	endif

c	Find the best FWHM
 11	fwhm = 2
	bestfwhm = fwhm
	minerror = 100000

 12	call smoothhist(n,n,hist2,fwhm,shist)
	call findlevels(n*n,shist,3,levels)
	call sensibleness(n,n,npts,hist,shist,levels,3,frac)
	frac(1) = (frac(1) - 68)**2 + (frac(2) - 95)**2 + (frac(3) -99)**2
	frac(1) = sqrt(frac(1))
	write(*,*) 'Find the best FWHM : ',fwhm,'/',n,' --> error = ',frac(1)
	if( frac(1).le.minerror ) then
	  minerror=frac(1)
	  bestfwhm=fwhm
	endif
	if( fwhm.lt.n ) then
	  fwhm = fwhm + 2
	  goto 12
	endif
	fwhm = bestfwhm
	write(*,*) 'Best FWHM : ',fwhm
	
 15	call smoothhist(n,n,hist2,fwhm,shist)
	call findlevels(n*n,shist,3,levels)
	call sensibleness(n,n,npts,hist,shist,levels,3,frac)
	smin = 0.0
	smax = -1e32
	do i=1,n
	  do j=1,n
	    smax = max(smax,shist(i,j))
	  enddo
	enddo
	
	write(*,'(a,f5.1,a)') 
     &'    68% contour contains ',frac(1),'% of samples,'
	write(*,'(a,f5.1,a)') 
     &'    95% contour contains ',frac(2),'% of samples,'
	write(*,'(a,f5.1,a)') 
     &'    99% contour contains ',frac(3),'% of samples,'
	write(*,*) 
	frac(1) = (frac(1) - 68)**2 + (frac(2) - 95)**2 + (frac(3) -99)**2
	frac(1) = sqrt(frac(1))
	write(*,'(a,f5.1)') 'error : ',frac(1)
	
c-- plot -------------------------------------------------------------

c Prepare the colortable and a square window
	call ctab (c,r,g,b,nc)
        if( plot.eq.1 ) then
		call pgbeg(0,'/xs',1,1)
		call pgask(.false.)
		call pgpage
		call pgpap(0.0,1.0)
		call pgscir(19,128)
		call pgctab(c,r,g,b,nc,1.,.5)
		call pgsch(1.2)
		call pgslw(3)
		
c		call pgqcir(c1,c2)
c		cdelt = 1./(c2 - c1)
c		do i=c1,(c2-c1)+c1
c		  call pgscr(i, 1.-cdelt*(i-c1)*cdelt*(i-c1),1.-cdelt*(i-c1),1.-sqrt(cdelt*(i-c1)))
c		enddo
c		call pgscr(0,1.,1.,1.)
c		call pgscr(1,0.,0.,0.)
c		do i=1,n
c		  do j=1,n
c		    shist(i,j)=( shist(i,j) - smin )/smax
c		    shist(i,j)=shist(i,j)*(c2-c1) + c1
c		  enddo
c		enddo
	
		call pgpage
		call pgvport(0.15,0.9,0.15,0.9)
		call pgwindow(xmin,xmax,ymin,ymax)
		call pgimag(shist,n,n,1,n,1,n,smin,smax,tr)
		call plotcontours(shist,n,tr,smin,smax)
		call pgsci(4)
		call pgpt(npts,x,y,1)
		call pgsci(1)
		call pgbox('BCTSN',0.0,0,'BCTSN',0.0,0)
		call pglab(label1,label2,title)
		
		write(*,'(a,$)') ' Re-smooth? (y or n): '
	       	read(*,'(a1)') ans
	        if (ans.eq.'y') then 
	       	  goto 10
	       	endif
	        write(8,'(a1)') 'n'
	       	write(8,*) fwhm 
	        write(8,'(a1)') 'n' 
	       	call pgend
        endif
c postscript output:

	write(*,*) 
	write(*,*) 'Plotting postscript...'
	write(*,*) 
	
	call pgbeg(0,'histogram2d.ps/vcps',1,1)
	call pgask(.false.)
	call pgpage
	call pgsclp(0)
	call pgpap(0.,1.)
	call pgsch(1.5)
	call pgslw(3)
	call pgscir(19,128)
	call pgctab(c,r,g,b,nc,1.,.5)
	call pgvport(0.15,0.9,0.2,0.9)
	call pgwindow(xmin,xmax,ymin,ymax)
	call pgimag(shist,n,n,1,n,1,n,smin,smax,tr)
	call plotcontours(shist,n,tr,smin,smax)
	write(*,'(a,$)') ' Plot points? (y or n): '
	read(*,'(a1)') ans
	write(8,'(a1)') ans 
	if (ans.eq.'y') then 
	  call pgsci(4)
c	  call pgsch(3.0)
	  if (npts.lt.100) then
            call pgpt(npts,x,y,17)
	  else
            call pgpt(npts,x,y,1)
	  endif
          call pgsci(1)
c	  call pgsch(1.5)
	endif
	write(*,*) 
	i = 1
 20	write(*,'(a,$)') ' Plot reference point? (y or n): '
c	Increment the color index i
	i = i + 1
	read(*,'(a1)') ans
	write(8,'(a1)') ans 
	if (ans.eq.'y') then
	  write(*,'(a,$)') '    Enter point (x,y): '
	  read(*,*) pt1,pt2
	  write(8,*) pt1,pt2 
	  write(*,'(a,$)') '    Enter errors on point (ex,ey): '
	  read(*,*) err1,err2
	  write(8,*) err1,err2
	  call pgscr(17,1.0,1.0,1.0)
	  call pgsci(i)
c	Change here the star size (originally 2.0)
	  call pgsch(4.0)
	  if (err1.eq.0.and.err2.eq.0) then
	    call pgsci(5.)
	    call pgpt1(pt1,pt2,18)
	  else 
	    call pgpt1(pt1,pt2,17)
	    call pgerr1(5,pt1,pt2,err1,1.0)
	    call pgerr1(6,pt1,pt2,err2,1.0)
	  endif
	  call pgsci(1)
	  call pgsch(1.2)
	  goto 20
	endif
	
	call pgsch(1.5)
	call pgbox('BCTSN',0.0,0,'BCTSN',0.0,0)
c	Change here the label size (originally 1.7)
	call pgsch(2.5)
	call pgmtxt('B',1.8,0.5,0.5,label1)
	call pgmtxt('L',1.3,0.5,0.5,label2)
	call pgmtxt('T',0.5,0.5,0.5,title)

	write(*,*) 
	call pgend
		
	end

c=======================================================================
	subroutine smoothhist(nx,ny,hist2,fwhm,shist)
	
	implicit none
	integer nx,ny
	real hist2(2*nx,2*ny),fwhm,shist(nx,ny)

	integer i,j,k,iw,nx2,ny2
	real shist2(2*nx,2*ny),sum

	nx2 = 2*nx
	ny2 = 2*ny

c	write(*,*) '  Smoothing histogram...'
c	write(*,*) 
	iw = 2
	do i=1,nx
	  do j=1,ny
	    shist(i,j) = 0.0
	  enddo
	enddo
	do i=1,nx2
	  do j=1,ny2
	    shist2(i,j) = 0.0
	  enddo
	enddo

	call smooth(hist2,nx2,ny2,shist2,nx2,ny2,nint(fwhm),iw)
	call crop(nx2,nx2,shist2,nx,ny,shist)

	sum = 0.0
	do i=1,nx
	  do j=1,ny
	    sum = sum + shist(i,j)
	  enddo
	enddo
	do i=1,nx
	  do j=1,ny
	    shist(i,j) = shist(i,j)/sum
	  enddo
	enddo
	
	return
	end

c==================================================================
c check sensibleness of contours
	subroutine sensibleness(nx,ny,npts,hist,shist,levels,nc, frac)

	implicit none
	integer nx,ny,npts,nc
	real hist(nx,ny),shist(nx,ny),levels(nc),frac(nc)

	integer i,j,k

	do k=1,nc
	  frac(k) = 0.0
	enddo
	do i=1,nx
	  do j=1,ny
	    do k=1,nc
	      if (shist(i,j).ge.levels(k)) then
	        frac(k) = frac(k) + hist(i,j)
	      endif
	    enddo
	  enddo
	enddo
	do k=1,nc
	  frac(k) = 100.0*frac(k)/float(npts)
	enddo
	
	return
	end
c=======================================================================
	
	subroutine pad(nx1,ny1,map1,nx2,ny2,map2)
	
	implicit none
	
	integer nx1,ny1,nx2,ny2
	real map1(nx1,ny1),map2(nx2,ny2)
	integer i,j,imargin,jmargin,ii,jj
	
      imargin = (nx2 - nx1 + 1)/2
      jmargin = (ny2 - ny1 + 1)/2

      do i=1,nx2
        do j=1,ny2
          map2(i,j) = 0.0
        enddo
      enddo
	
      do i=1,nx1
        do j=1,ny1
          ii = i + imargin 
          jj = j + jmargin 
          map2(ii,jj) = map1(i,j)
        enddo
      enddo
	
	return
	end
	
c=======================================================================
	
	subroutine crop(nx2,ny2,map2,nx1,ny1,map1)
	
	implicit none
	
	integer nx1,ny1,nx2,ny2
	real map1(nx1,ny1),map2(nx2,ny2)
	integer i,j,imargin,jmargin,ii,jj
	
      imargin = (nx2 - nx1 + 1)/2
      jmargin = (ny2 - ny1 + 1)/2

      do i=1,nx1
        do j=1,ny1
          ii = i + imargin 
          jj = j + jmargin 
          map1(i,j) = map2(ii,jj)
        enddo
      enddo
            
	return
	end	
	
c=======================================================================


      subroutine ctab(x,r,g,b,nc)

      integer nc
      real x(nc),r(nc),g(nc),b(nc)

      x(1) = 0.0
      x(2) = 0.25
      x(3) = 0.5
      x(4) = 0.75
      x(5) = 1.0

c Colours set as linear greyscale 

      r(1)=1.00
      g(1)=1.00
      b(1)=1.00

      r(2)=0.85
      g(2)=0.85
      b(2)=0.85

      r(3)=0.7
      g(3)=0.7
      b(3)=0.7

      r(4)=0.55
      g(4)=0.55
      b(4)=0.55

      r(5)=0.4
      g(5)=0.4
      b(5)=0.4

cc Colours set as black body orange spectrum
cc
cc      r(5)=0.00
cc      g(5)=0.00
cc      b(5)=0.00
cc
cc      r(4)=0.5
cc      g(4)=0.0
cc      b(4)=0.0
cc
cc      r(3)=1.0
cc      g(3)=0.5
cc      b(3)=0.0
cc
cc      r(2)=1.0
cc      g(2)=1.0
cc      b(2)=0.5
cc
cc      r(1)=1.0
cc      g(1)=1.0
cc      b(1)=1.0

	return
	end
	
c=======================================================================

      subroutine plotcontours(dummy,size,trans,min,max)

      implicit none

      integer size
      real dummy(size,size),trans(6),min,max
      integer nc
      parameter (nc=3)
      real c(nc)
       
c      do i=1,nc
c        c(i) = min + (max-min)*float(i-1)/float(nc-1)
c      enddo
  
	if (min.eq.0.0) then 
	  call findlevels(size*size,dummy,nc,c)
	else
	  call finddblevels(size*size,dummy,nc,c)
	endif
	
      call pgscr(17,1.0,1.0,1.0)
      call pgscr(18,0.0,0.0,0.0)
      call pgslw(3)
c      call pgsci(17)
c      call pgcont(dummy,size,size,1,size,1,size,c,-nc,trans)                 
c      call pgslw(1)
      call pgsci(18)
      call pgcont(dummy,size,size,1,size,1,size,c,-nc,trans)     
      call pgsci(1)
c       call pgslw(1)
      
	return
      end
c
c=======================================================================       
c
        subroutine bin(nx,ny,m,x,y,z,xcell,ycell,tr)
c
        implicit none
c       
        integer nx,ny,m
        real x(m),y(m),z(nx,ny)
c       
        real xcell,ycell,tr(6)
c
        real xmin,xmax,ymin,ymax,xa,xb,ya,yb,sum      
      integer count
      integer i,j,icentre,jcentre,k
c
c-----------------------------------------------------------------------

        sum = 0.0
	  do i=1,nx
            do j=1,ny
c
              xa = tr(1) + tr(2)*i - xcell/2.0
              xb = xa + xcell    
              ya = tr(4) + tr(6)*j + ycell/2.0
              yb = ya + ycell    
c         
              count = 0
              do k=1,m
                if (x(k).ge.xa.and.x(k).lt.xb
     &                  .and.y(k).ge.ya.and.y(k).lt.yb) then
                    count = count + 1
                endif
              enddo

              z(i,j) = float(count)
           
            enddo
          enddo
        
	  
	return
	end

c End of binning shenanigans    
c-----------------------------------------------------------------------

	include 'smooth.f'
	
c-----------------------------------------------------------------------

	subroutine findlevels(n,z,nc,levels)
	
	integer n,nc
	real z(n),copy(n),levels(nc)
	
	real zmax,sum
	integer i,count,current
	
	sum = 0.0
	do i=1,n
	  copy(i) = z(i)
	  sum = sum + copy(i)
	enddo
	
	call sort(copy,n)
	
	current = 1
	sum = 0.0
	do i=1,n
	  sum = sum + copy(i)
	  if (sum.ge.0.68.and.current.eq.1.and.nc.ge.1) then
	    levels(1) = copy(i)
	    current = 2
	  elseif (sum.ge.0.95.and.current.eq.2.and.nc.ge.2) then
	    levels(2) = copy(i)
	    current = 3
	  elseif (sum.ge.0.99.and.current.eq.3.and.nc.ge.3) then
	    levels(3) = copy(i)
	    current = 4
	  endif
	  if (current.eq.4) goto 999
	enddo
	
	
	
 999	return
	end
	
c-----------------------------------------------------------------------

	subroutine finddblevels(n,z,nc,levels)
	
	integer n,nc
	real z(n),copy(n),levels(nc)
	
	real zmax,sum
	integer i,count,current
	
	sum = 0.0
	do i=1,n
	  copy(i) = 10.0**(z(i)/10.0)
	  sum = sum + copy(i)
	enddo
	
	call sort(copy,n)
	
	current = 1
	sum = 0.0
	do i=1,n
	  sum = sum + copy(i)
	  if (sum.ge.0.68.and.current.eq.1.and.nc.ge.1) then
	    levels(1) = copy(i)
	    current = 2
	  elseif (sum.ge.0.95.and.current.eq.2.and.nc.ge.2) then
	    levels(2) = copy(i)
	    current = 3
	  elseif (sum.ge.0.99.and.current.eq.3.and.nc.ge.3) then
	    levels(3) = copy(i)
	    current = 4
	  endif
	  if (current.eq.4) goto 999
	enddo
	
 999	do i=1,nc
	  levels(i) = 10*log10(levels(i)) 
	enddo
	
	return
	end

!=======================================================================

	include 'sort.f'

!=======================================================================