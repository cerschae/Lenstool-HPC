c-----*-----------------------------------------------------------------

      subroutine regridarray(nx,ny,nz,ninx,niny,noutx,nouty,array,
     &                       dummy,dummy2,xp,yp,iway)
      implicit none
      integer  nx,ny,nz,iway
      integer  ninx(nz),niny(nz),noutx(nz),nouty(nz)
      real     array(nx,ny,nz)
      real     dummy(nx,ny),dummy2(nx,ny),xp(nx),yp(ny)

      integer  i,j,m,nxin,nyin,nxout,nyout

c ... begin loop over last index
      do m=1,nz
        nxin=ninx(m)
        nyin=niny(m)
        nxout=noutx(m)
        nyout=nouty(m)
c ... copy relevant data into dummy
        do i=1,nx
          do j=1,ny
            if (i.le.nxin .and. j.le.nyin) then
              dummy(i,j)=array(i,j,m)
            else
              dummy(i,j)=0.
            end if
          end do
        end do
c ... regrid dummy
        call regrid(nx,ny,nxin,nyin,nxout,nyout,dummy,dummy2,xp,yp,iway)
c ... copy back into array 
        do i=1,nx
          do j=1,ny
            array(i,j,m)=dummy2(i,j)
          end do
        end do
c ... close loop over last index
      end do

      return
      end

c-----*-----------------------------------------------------------------

c subroutine to regrid array using interp

      subroutine regrid(nx,ny,nxin,nyin,nxout,nyout,arrin,arrout,xp,yp,
     &                  iw)
      implicit none
      integer  nx,ny,nxin,nyin,nxout,nyout,iw
      real     arrin(nx,ny),arrout(nx,ny)
      real     xp(nx),yp(ny)

      integer  i,j
      real     xnew,ynew,value

c ... if input and output sizes equal then copy real part and return
      if (nxin.eq.nxout .and. nyin.eq.nyout) then
        do i=1,nx
          do j=1,ny
            arrout(i,j)=arrin(i,j)
          end do
        end do
        return
      end if
c ... create arrays of measured x and y values
      do i=1,nxin
        xp(i)=0.5 + (float(2*i-1)/2.0)*float(nx)/float(nxin)
      end do
      do i=1,nyin
        yp(i)=0.5 + (float(2*i-1)/2.0)*float(ny)/float(nyin)
      end do
c ... interpolate array
      do i=1,nx
        do j=1,ny
          if (i.le.nxout .and. j.le.nyout) then
            xnew=0.5 + (float(2*i-1)/2.0)*float(nx)/float(nxout)
            ynew=0.5 + (float(2*j-1)/2.0)*float(ny)/float(nyout)
            call interp(nx,ny,nxin,nyin,arrin,xp,yp,xnew,ynew,value,iw)
            arrout(i,j)=value
          else
            arrout(i,j)=0.
          end if
        end do
      end do

      end

c-----*-----------------------------------------------------------------

c subroutine to perform bilinear or nearest-neighbour interpolation

      subroutine interp(nx,ny,nxin,nyin,image,xp,yp,xnew,ynew,value,iw)
      implicit none
      integer  nx,ny,nxin,nyin,iw
      real     image(nx,ny)
      real     xp(nx),yp(ny)
      real     xnew,ynew,value

      integer  i,j,flag
      real     xoldstep,yoldstep,t,u,bl,br,tr,tl

      flag=0
      xoldstep=xp(2)-xp(1)
      yoldstep=yp(2)-yp(1)
      i=int((xnew-0.5)*float(nxin)/float(nx)+0.5)
      j=int((ynew-0.5)*float(nyin)/float(ny)+0.5)
      t=(xnew-xp(i))/xoldstep
      u=(ynew-yp(j))/yoldstep

      if (i.lt.1) then
        i=1
        flag=1
      end if
      if (i.ge.nxin) then
        i=nxin
        flag=1
      end if
      if (j.lt.1) then
        j=1
        flag=1
      end if
      if (j.ge.nyin) then
        j=nyin
        flag=1
      end if

      if (iw.eq.0) then
c ... perform bilinear interpolation
        if (flag.eq.1) then
          value=image(i,j)
        else
          bl=image(i,j)
          br=image(i+1,j)
          tr=image(i+1,j+1)
          tl=image(i,j+1)
          value=(1.0-t)*(1.0-u)*bl
     &         +t*(1.0-u)*br
     &         +t*u*tr
     &         +(1.0-t)*u*tl
        end if
      else if (iw.eq.1) then
c ... perform nearest neighbour interpolation
        if (flag.eq.1) then
          value=image(i,j)
        else
          bl=image(i,j)
          br=image(i+1,j)
          tr=image(i+1,j+1)
          tl=image(i,j+1)
          if (t.lt.0.5) then
            if (u.lt.0.5) value=bl
            if (u.ge.0.5) value=tl
          else
            if (u.lt.0.5) value=br
            if (u.ge.0.5) value=tr
          end if
        end if
      end if

      return
      end

c-----*-----------------------------------------------------------------
