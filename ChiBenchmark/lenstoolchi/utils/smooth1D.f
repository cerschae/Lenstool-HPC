c=======================================================================
c
      subroutine smooth1D(small,nxin,array,nx,npix,iw)
c
c Convolves an nxin array with a gaussian of width bfwhm pixels,
c writing the result to a nx array. 
c   iw = 0  ->  bilinear interpolation (iw = 2 -> silent)
c   iw = 1  ->  nearest neighbours interp. (iw = 3 -> silent)
c
      implicit none
c	
      integer nxin,nx,nxout
	integer npix,iw
	real small(nxin),large(nx)
      real array(nx)
      real value
      real barea,bfwhm
      real beam(nx),work(2*nx)
      real xp(nx)
      complex beamfft(nx),wkspce(nx),wkspce2(nx)
	logical verbose
c      
      integer i,j,ndim,nn(2),n
      integer job,iform
c
c-----------------------------------------------------------------------
c
c initialise variables

	nxout = nx
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
	  large(i) = 0.0
	enddo    
      do i=1,nxin
	  large(i) = small(i)
	enddo    	
	
	ndim=1
      nn(1)=nx

      if (verbose) write(*,*) 'Regridding array ... '
      call regrid1D(nx,nxin,nxout,large,array,xp,iw)

c ... make beam

      if (verbose) write(*,*) 'Making Gaussian ... '
      call makebeam(beam,nx,nxout,barea,bfwhm)

c ... make beam FFT

      n=0
      do i=1,nxout
        n=n+1
        wkspce(n)=cmplx(beam(i),0.)
      end do
      job=-1                   ! forward transform (-ve exponential)
      iform=0                  ! data are real
      call cheqboard(wkspce,nxout)
      call fourt(wkspce,nn,ndim,job,iform,work)
      call cheqboard(wkspce,nxout)
      n=0
      do i=1,nxout
        n=n+1
        beamfft(i)=wkspce(n)
      end do

c ... convolve maps with beam

      if (verbose) write(*,*) 'Convolving map with Gaussian ...'
      n=0
      do i=1,nxout
        n=n+1
        value=array(i)
        wkspce(n)=cmplx(value,0.)
      end do
	
c .. FFT sky
      job=-1                     ! forward transform (-ve exponential)
      iform=0                    ! data are real
	
      call cheqboard(wkspce,nxout)
      call fourt(wkspce,nn,ndim,job,iform,work)
      call cheqboard(wkspce,nxout)

c ... multiply FFT(sky) by beamfft at given frequency

      n=0
      do i=1,nxout
        n=n+1
        wkspce2(n)=wkspce(n)*beamfft(i)
      end do
	
c ... inverse FFT to obtain convolved sky at given frequency

      job=1                      ! backward transform (+ve exponential)
      iform=1                    ! data are complex
      call cheqboard(wkspce2,nxout)
      call fourt(wkspce2,nn,ndim,job,iform,work)
      call cheqboard(wkspce2,nxout)
	
c ... normalise FFT and copy to array

      n=0
      do i=1,nxout
        n=n+1
        array(i)=real(wkspce2(n))/real(nx)/barea
      end do

 999	return
      end

c-----*-----------------------------------------------------------------
 
      subroutine cheqboard(array,nx)
      implicit none
 
      integer nx
      complex array(nx)
 
      integer i,j
 
      do i=1,nx-1,2
        array(i+1)=-array(i+1)
      end do
 
      end
 
c-----*-----------------------------------------------------------------

      subroutine makebeam(beam,nx,nxout,beamarea,fwhm)

      implicit none
      integer  nx,nxout
      real     beam(nx),beamarea,fwhm

      integer  i,j
      real     sigma,dist
      integer  centrex

      do i=1,nx
        beam(i)=0.
      end do

      centrex = nxout/2 + 1
      sigma=fwhm/(2.*sqrt(2.*log(2.)))
 
      beamarea=0.
      do i=1,nxout
        dist=sqrt(float(i-centrex)**2)
        beam(i)=exp(-dist**2 / (2*sigma**2))
        beamarea=beamarea+beam(i)
      end do
 
      return
      end

c-----*-----------------------------------------------------------------

      include "rfft.f"
	include "regrid1D.f"
	
c-----*-----------------------------------------------------------------


