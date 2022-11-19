!*******************************************************************************
! this can have those function to generate evenly spaced and non-uniformly spaced nodes, 
! like circular distributed nodes and arbitrarily distributed.

subroutine generate1DPts(x0,x1,n,xpts)
! this is for generating evenly spaced points
	implicit none

	real,intent(in)::x0,x1
	integer,intent(in)::n
	real,intent(out)::xpts(n)

	integer::i
	real::dx
	
	! Initialization
	xpts = dble(0)

	dx = (x1-x0)/(n-1)

	do i=1,n
		xpts(i) = x0 + dx*(i-1)
	end do

end subroutine

!*******************************************************************************

!*******************************************************************************

subroutine generate2dPts(x_start,x_end,npx,y_start,y_end,npy,xynodes)
	implicit none
	
	real,intent(in)::x_start,x_end,y_start,y_end
	integer,intent(in)::npx,npy
	real,intent(out)::xynodes(npx*npy,2)
	
	real::xpts(npx),ypts(npy)
	integer::i,j,k
	
	write(*,*)"Generating the nodal points..."
	
	xynodes = dble(0)
	
	call generate1DPts(x_start,x_end,npx,xpts)
	call generate1DPts(y_start,y_end,npy,ypts)
	
	do i = 1,npy
		do j = 1,npx
			k = (i-1)*npx+j
			xynodes(k,1) = xpts(j)
			xynodes(k,2) = ypts(i)
		end do 
	end do 
	
	write(*,*)"Done!"
	
end subroutine
	
!*******************************************************************************