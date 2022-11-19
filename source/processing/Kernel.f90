!*******************************************************************************

subroutine Kernel1D(phi, phidx, xpi, xx, h)
! 	This is to calculate the kernel function and its first derivative
!   phi(x) = phi((xx-xpi)/h)
! 	using cubic spline kernel function
	implicit none

	real,intent(in)::xx,xpi,h
	real,intent(out)::phi,phidx
	
	real::z
	z = (xx-xpi)/h
	
	phi = dble(0)
	phidx = dble(0)

	if (-1. <= z .and. z <= -0.5) then
		phi = 4./3. + 4.*z + 4.*z**2 + 4./3.*z**3
		phidx = 4.*(1./h) + 8.*z*(1./h) + 4.*z**2*(1./h)

	elseif (-.5 <= z .and. z <= 0.) then
		phi = 2./3. - 4.*z**2 - 4.*z**3
		phidx = -8.*z*(1./h) - 12.*z**2*(1/h)

	elseif ( 0 <= z .and. z <= .5) then
		phi = 2./3. - 4.*z**2 + 4.*z**3
		phidx = -8.*z*(1./h) + 12.*z**2*(1/h)

	elseif( 0.5 <= z .and. z <= 1. ) then
		phi = 4./3. - 4.*z + 4.*z**2 - 4./3.*z**3
		phidx = -4.*(1./h) + 8.*z*(1./h) - 4.*z**2*(1./h)

	else
		phi = 0.
		phidx = 0.
	end if

end subroutine

!*******************************************************************************

subroutine Kernel2D_rec(phi, phidxy, xpi, ypi, xx, yy)
	use constants,only:hx,hy
	implicit none
	
	real,intent(in)::xx,yy,xpi,ypi
	real,intent(out)::phi,phidxy(2) ! (fx,fy)
	
	real::phix, phidx, phiy, phidy
	
	phi = dble(0)
	phidxy = (/ dble(0), dble(0) /)
	phix = dble(0)
	phidx = dble(0)
	phiy = dble(0)
	phidy = dble(0)
	
	call Kernel1D(phix, phidx, xpi, xx, hx)
	call Kernel1D(phiy, phidy, ypi, yy, hy)
	
	phi = phix*phiy 
	phidxy(1) = phidx*phiy
	phidxy(2) = phix*phidy

end subroutine

!*******************************************************************************

subroutine Kernel2D_cir(phi, phidxy, xpi, ypi, xx, yy)
	use constants
	implicit none
	
	real,intent(in)::xx,yy,xpi,ypi
	real,intent(out)::phi,phidxy(2) ! (fx,fy)
	
	real::z,phidz
	
	phi = dble(0)
	phidz = dble(0)
	phidxy = (/ dble(0), dble(0) /)
	
	z = sqrt((xx-xpi)**2 + (yy-ypi)**2)/hr
	
	if (-1. <= z .and. z <= -0.5) then
		phi = 4./3. + 4.*z + 4.*z**2 + 4./3.*z**3
		phidz = 4. + 8.*z + 4.*z**2

	elseif (-.5 <= z .and. z <= 0.) then
		phi = 2./3. - 4.*z**2 - 4.*z**3
		phidz = -8.*z - 12.*z**2

	elseif ( 0 <= z .and. z <= .5) then
		phi = 2./3. - 4.*z**2 + 4.*z**3
		phidz = -8.*z + 12.*z**2

	elseif( 0.5 <= z .and. z <= 1. ) then
		phi = 4./3. - 4.*z + 4.*z**2 - 4./3.*z**3
		phidz = -4. + 8.*z - 4.*z**2

	else
		phi = 0.
		phidz = 0.
	end if
	
	phidxy(1) = phidz*(xx-xpi)/(hr*z)
	phidxy(2) = phidz*(yy-ypi)/(hr*z)	

end subroutine

!*******************************************************************************