c=======================================================================
c
      subroutine smooth(small,nxin,nyin,array,nx,ny,npix,iw)
c
c Convolves an nxin by nyin image with a gaussian of width npix pixels,
c writing the result to a nx by ny array. 
c   iw = 0  ->  bilinear interpolation (iw = 2 -> silent)
c   iw = 1  ->  nearest neighbours interp. (iw = 3 -> silent)
c
      implicit none
c	
      integer nxin,nyin,nx,ny,nxout,nyout
	integer npix,iw
	real small(nxin,nyin),large(nx,ny)
      real array(nx,ny)
      real value
      real barea,bfwhm
      real beam(nx,ny),work(2*nx)
      real xp(nx),yp(ny)
      complex beamfft(nx,ny),wkspce(nx*ny),wkspce2(nx*ny)
	logical verbose
c      
      integer i,j,ndim,nn(2),n
      integer job,iform
c
c-----------------------------------------------------------------------
c
c initialise variables

	nxout = nx
	nyout = ny
      bfwhm = float(npix)
	if (iw.eq.2) then
	  iw = 0
	  verbose = .false.
	elseif (iw.eq.3) then  
	  iw = 1
	  verbose = .false.
	else 
	  verbose = .true.
	endif    
	
c regrid array

      do i=1,nx
	  do j=1,ny
	    large(i,j) = 0.0
	  enddo
	enddo    
      do i=1,nxin
	  do j=1,nyin
	    large(i,j) = small(i,j)
	  enddo
	enddo    	
	
	ndim=2
      nn(1)=nx
      nn(2)=ny

      if (verbose) write(*,*) 'Regridding array ... '
      call regrid(nx,ny,nxin,nyin,nxout,nyout,large,array,xp,yp,iw)

c ... make beam

      if (verbose) write(*,*) 'Making Gaussian ... '
      call makebeam(beam,nx,ny,nxout,nyout,barea,bfwhm)

c ... make beam FFT

      n=0
      do i=1,nxout
        do j=1,nyout
          n=n+1
          wkspce(n)=cmplx(beam(i,j),0.)
        end do
      end do
      job=-1                   ! forward transform (-ve exponential)
      iform=0                  ! data are real
      call cheqboard(wkspce,nxout,nyout)
      call fourt(wkspce,nn,ndim,job,iform,work)
      call cheqboard(wkspce,nxout,nyout)
      n=0
      do i=1,nxout
        do j=1,nyout
          n=n+1
          beamfft(i,j)=wkspce(n)
        end do
      end do

c ... convolve maps with beam

      if (verbose) write(*,*) 'Convolving map with Gaussian ...'
      n=0
      do i=1,nxout
        do j=1,nyout
          n=n+1
          value=array(i,j)
          wkspce(n)=cmplx(value,0.)
        end do
      end do
	
c .. FFT sky
      job=-1                     ! forward transform (-ve exponential)
      iform=0                    ! data are real
	
      call cheqboard(wkspce,nxout,nyout)
      call fourt(wkspce,nn,ndim,job,iform,work)
      call cheqboard(wkspce,nxout,nyout)

c ... multiply FFT(sky) by beamfft at given frequency

      n=0
      do i=1,nxout
        do j=1,nyout
          n=n+1
          wkspce2(n)=wkspce(n)*beamfft(i,j)
        end do
      end do
	
c ... inverse FFT to obtain convolved sky at given frequency

      job=1                      ! backward transform (+ve exponential)
      iform=1                    ! data are complex
      call cheqboard(wkspce2,nxout,nyout)
      call fourt(wkspce2,nn,ndim,job,iform,work)
      call cheqboard(wkspce2,nxout,nyout)
	
c ... normalise FFT and copy to array

      n=0
      do i=1,nxout
        do j=1,nyout
          n=n+1
          array(i,j)=real(wkspce2(n))/real(nx*ny)/barea
        end do
      end do

 999	return
      end

c-----*-----------------------------------------------------------------
 
      subroutine cheqboard(array,nx,ny)
      implicit none
 
      integer nx,ny
      complex array(nx,ny)
 
      integer i,j
 
      do i=1,nx-1,2
        do j=1,ny-1,2
          array(i+1,j)=-array(i+1,j)
          array(i,j+1)=-array(i,j+1)
        end do
      end do
 
      end
 
c-----*-----------------------------------------------------------------

      subroutine makebeam(beam,nx,ny,nxout,nyout,beamarea,fwhm)

      implicit none
      integer  nx,ny,nxout,nyout
      real     beam(nx,ny),beamarea,fwhm

      integer  i,j
      real     sigma,dist
      integer  centrex,centrey

      do i=1,nx
        do j=1,ny
          beam(i,j)=0.
        end do
      end do

      centrex = nxout/2 + 1
      centrey = nyout/2 + 1
      sigma=fwhm/(2.*sqrt(2.*log(2.)))
 
      beamarea=0.
      do i=1,nxout
        do j=1,nyout
          dist=sqrt(float(i-centrex)**2+float(j-centrey)**2)
          beam(i,j)=exp(-dist**2 / (2*sigma**2))
          beamarea=beamarea+beam(i,j)
        end do
      end do
 
      return
      end

c-----*-----------------------------------------------------------------

      include "rfft.f"
	include "regrid.f"
	
c-----*-----------------------------------------------------------------


