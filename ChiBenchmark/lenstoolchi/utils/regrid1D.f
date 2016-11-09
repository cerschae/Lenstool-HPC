cc-----*-----------------------------------------------------------------

c subroutine to regrid array using interp

      subroutine regrid1D(nx,nxin,nxout,arrin,arrout,xp,iw)
      implicit none
      integer  nx,nxin,nxout,iw
      real     arrin(nx),arrout(nx)
      real     xp(nx)

      integer  i,j
      real     xnew,value

c ... if input and output sizes equal then copy real part and return
      if (nxin.eq.nxout) then
        do i=1,nx
          arrout(i)=arrin(i)
        end do
        return
      end if
c ... create arrays of measured x and y values
      do i=1,nxin
        xp(i)=0.5 + (float(2*i-1)/2.0)*float(nx)/float(nxin)
      end do
c ... interpolate array
      do i=1,nx
        if (i.le.nxout) then
          xnew=0.5 + (float(2*i-1)/2.0)*float(nx)/float(nxout)
          call interp1D(nx,nxin,arrin,xp,xnew,value,iw)
          arrout(i)=value
        else
          arrout(i)=0.
        end if
      end do

      end

c-----*-----------------------------------------------------------------

c subroutine to perform bilinear or nearest-neighbour interpolation

      subroutine interp1D(nx,nxin,image,xp,xnew,value,iw)
      implicit none
      integer  nx,nxin,iw
      real     image(nx)
      real     xp(nx)
      real     xnew,value

      integer  i,j,flag
      real     xoldstep,t,u,bl,br,tr,tl

      flag=0
      xoldstep=xp(2)-xp(1)
      i=int((xnew-0.5)*float(nxin)/float(nx)+0.5)
      t=(xnew-xp(i))/xoldstep

      if (i.lt.1) then
        i=1
        flag=1
      end if
      if (i.ge.nxin) then
        i=nxin
        flag=1
      end if

      if (iw.eq.0) then
c ... perform linear interpolation
        if (flag.eq.1) then
          value=image(i)
        else
          bl=image(i)
          br=image(i+1)
          value=(1.0-t)*bl + t*br
        end if
      else if (iw.eq.1) then
c ... perform nearest neighbour interpolation
        if (flag.eq.1) then
          value=image(i)
        else
          bl=image(i)
          br=image(i+1)
          if (t.lt.0.5) then
            value=bl
          else
            value=br
          end if
        end if
      end if

      return
      end

c-----*-----------------------------------------------------------------
