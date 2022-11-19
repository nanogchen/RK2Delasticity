!*******************************************************************************
! This should include multi-base of P. And use the matrix operation to calculate M(x),N(x)

subroutine CalN(N, Ndxy, xynodes, xpi, ypi, xx, yy)
	use constants
	implicit none
	
	real,intent(in)::xpi, ypi, xx, yy
	real,intent(in)::xynodes(np,2)
	real,intent(out)::N, Ndxy(2)
	
	real::M(3,3),Mdx(3,3),Mdy(3,3),invM(3,3),invMdx(3,3),invMdy(3,3) ! moment matrix and its inverse, derivative
	real::P(3), Pdx(3), Pdy(3),PT0(3)
	real::phi, phidxy(2)
	real::m1,m2,m3,m4,m5,m6
	real::m1dxy(2),m2dxy(2),m3dxy(2),m4dxy(2),m5dxy(2),m6dxy(2)
	
	real::xpii,ypii,Ndxy_temp1,Ndxy_temp2,Ndxy_temp3
	real::Ptemp(3),Ptemp1(3),Ptemp2(3),Ptemp3(3)
	integer::i,j,k

	! initialization	
	N = dble(0)
	Ndxy = (/ dble(0), dble(0) /)
	
	phi = dble(0)
	phidxy = (/ dble(0),dble(0)/)
	
	Ndxy_temp1 = dble(0)
	Ndxy_temp2 = dble(0)
	Ndxy_temp3 = dble(0)
	
	PT0 = (/ dble(1), dble(0), dble(0) /)
	Ptemp = (/ dble(0), dble(0), dble(0) /)
	Ptemp1 = (/ dble(0), dble(0), dble(0) /)
	Ptemp2 = (/ dble(0), dble(0), dble(0) /)
	Ptemp3 = (/ dble(0), dble(0), dble(0) /)
	P(1) = dble(1)
	P(2) = (xx - xpi)/hx
	P(3) = (yy - ypi)/hy
	Pdx(1) = dble(0)
	Pdx(2) = 1./hx
	Pdx(3) = dble(0)
	Pdy(1) = dble(0)
	Pdy(2) = dble(0)
	Pdy(3) = 1./hy
	
	m1 = dble(0)
	m2 = dble(0)
	m3 = dble(0)
	m4 = dble(0)
	m5 = dble(0)
	m6 = dble(0)
	m1dxy = (/dble(0),dble(0)/)
	m2dxy = (/dble(0),dble(0)/)
	m3dxy = (/dble(0),dble(0)/)
	m4dxy = (/dble(0),dble(0)/)
	m5dxy = (/dble(0),dble(0)/)
	m6dxy = (/dble(0),dble(0)/)
	
	do i = 1,3
		do j = 1,3
			M(i,j) = dble(0)			
			Mdx(i,j) = dble(0)
			Mdy(i,j) = dble(0)
			invM(i,j) = dble(0)
			invMdx(i,j) = dble(0)
			invMdy(i,j) = dble(0)
		end do
	end do
	
	do i=1,npy
		do j = 1,npx
			k = (i-1)*npx + j
			xpii = xynodes(k,1)
			ypii = xynodes(k,2)
			call Kernel2D_rec(phi, phidxy, xpii, ypii, xx, yy)
			m1 = m1 + phi
			m2 = m2 + phi*((xx-xpii)/hx)**2
			m3 = m3 + phi*((yy-ypii)/hy)**2
			m4 = m4 + phi*((xx-xpii)/hx)*((yy-ypii)/hy)
			m5 = m5 + phi*((yy-ypii)/hy)
			m6 = m6 + phi*((xx-xpii)/hx)
			
			m1dxy(1) = m1dxy(1) + phidxy(1)
			m2dxy(1) = m2dxy(1) + 2*((xx-xpii)/hx)*(1./hx)*phi + ((xx-xpii)/hx)**2*phidxy(1)
			m3dxy(1) = m3dxy(1) + ((yy-ypii)/hy)**2*phidxy(1)
			m4dxy(1) = m4dxy(1) + ((yy-ypii)/hy)*(1./hx)*phi + ((xx-xpii)/hx)*((yy-ypii)/hy)*phidxy(1)
			m5dxy(1) = m5dxy(1) + ((yy-ypii)/hy)*phidxy(1)
			m6dxy(1) = m6dxy(1) + (1./hx)*phi + ((xx-xpii)/hx)*phidxy(1)
			
			m1dxy(2) = m1dxy(2) + phidxy(2)
			m2dxy(2) = m2dxy(2) + ((xx-xpii)/hx)**2*phidxy(2)
			m3dxy(2) = m3dxy(2) + 2*((yy-ypii)/hy)*(1./hy)*phi + ((yy-ypii)/hy)**2*phidxy(2)
			m4dxy(2) = m4dxy(2) + ((xx-xpii)/hx)*(1./hy)*phi + ((xx-xpii)/hx)*((yy-ypii)/hy)*phidxy(2)
			m5dxy(2) = m5dxy(2) + (1./hy)*phi + ((yy-ypii)/hy)*phidxy(2)
			m6dxy(2) = m6dxy(2) + ((xx-xpii)/hx)*phidxy(2)		
		end do
	end do
	
	M(1,1) = m1
	M(1,2) = m6
	M(1,3) = m5
	M(2,1) = m6
	M(2,2) = m2
	M(2,3) = m4
	M(3,1) = m5
	M(3,2) = m4
	M(3,3) = m3
	
	call INVERSE(M,3,invM) 
		
	Mdx(1,1) = m1dxy(1)
	Mdx(1,2) = m6dxy(1)
	Mdx(1,3) = m5dxy(1)
	Mdx(2,1) = m6dxy(1)
	Mdx(2,2) = m2dxy(1)
	Mdx(2,3) = m4dxy(1)
	Mdx(3,1) = m5dxy(1)
	Mdx(3,2) = m4dxy(1)
	Mdx(3,3) = m3dxy(1)
	
	Mdy(1,1) = m1dxy(2)
	Mdy(1,2) = m6dxy(2)
	Mdy(1,3) = m5dxy(2)
	Mdy(2,1) = m6dxy(2)
	Mdy(2,2) = m2dxy(2)
	Mdy(2,3) = m4dxy(2)
	Mdy(3,1) = m5dxy(2)
	Mdy(3,2) = m4dxy(2)
	Mdy(3,3) = m3dxy(2)
	
	invMdx = matmul(invM,Mdx)
	invMdx = matmul(invMdx,invM)
	invMdx = -1.*invMdx
	
	invMdy = matmul(invM,Mdy)
	invMdy = matmul(invMdy,invM)
	invMdy = -1.*invMdy
	
	call Kernel2D_rec(phi, phidxy, xpi, ypi, xx, yy)
	
	Ptemp = matmul(PT0,invM)
	do i=1,3
		N = N + Ptemp(i)*P(i)
	end do
	! N = matmul(Ptemp,P)
	!******************************************************************************
	N = N*phi ! shape function
	
	Ptemp1 = matmul(PT0,invMdx)
	do i = 1,3
		Ndxy_temp1 = Ndxy_temp1 + Ptemp1(i)*P(i)
	end do
	! Ptemp1 = matmul(PT0,invMdy)
	Ndxy_temp1 = Ndxy_temp1*phi
	
	Ptemp2 = matmul(PT0,invM)
	do i = 1,3
		Ndxy_temp2 = Ndxy_temp2 + Ptemp2(i)*Pdx(i)
	end do
	! Ndxy_temp2 = matmul(Ptemp2,Pdx)
	Ndxy_temp2 = Ndxy_temp2*phi
	
	Ptemp3 = matmul(PT0,invM)
	do i = 1,3
		Ndxy_temp3 = Ndxy_temp3 + Ptemp3(i)*P(i)
	end do
	! Ndxy_temp3 = matmul(Ptemp3,P)
	Ndxy_temp3 = Ndxy_temp3*phidxy(1)
	
	Ndxy(1) = Ndxy_temp1 + Ndxy_temp2 + Ndxy_temp3 ! derivative of shape function
	!******************************************************************************
	
	!reinitialization
	Ndxy_temp1 = dble(0)
	Ndxy_temp2 = dble(0)
	Ndxy_temp3 = dble(0)
	
	Ptemp1 = (/ dble(0), dble(0), dble(0) /)
	Ptemp2 = (/ dble(0), dble(0), dble(0) /)
	Ptemp3 = (/ dble(0), dble(0), dble(0) /)
	
	!******************************************************************************
	Ptemp1 = matmul(PT0,invMdy)
	do i = 1,3
		Ndxy_temp1 = Ndxy_temp1 + Ptemp1(i)*P(i)
	end do
	! Ndxy_temp1 = matmul(Ptemp1,P)
	Ndxy_temp1 = Ndxy_temp1*phi
	
	Ptemp2 = matmul(PT0,invM)
	do i = 1,3
		Ndxy_temp2 = Ndxy_temp2 + Ptemp2(i)*Pdy(i)
	end do
	! Ndxy_temp2 = matmul(Ptemp2,Pdy)
	Ndxy_temp2 = Ndxy_temp2*phi
	
	Ptemp3 = matmul(PT0,invM)
	do i = 1,3
		Ndxy_temp3 = Ndxy_temp3 + Ptemp3(i)*P(i)
	end do
	! Ndxy_temp3 = matmul(Ptemp3,P)
	Ndxy_temp3 = Ndxy_temp3*phidxy(2)
	
	Ndxy(2) = Ndxy_temp1 + Ndxy_temp2 + Ndxy_temp3 ! derivative of shape function
	!******************************************************************************
	
end subroutine

!*******************************************************************************

subroutine CalN_FE(Nfe,pti,pt)
	use constants
!	This is to construct the finite element shape function along the essential boundary
	implicit none
	
	real,intent(in)::pti(2),pt(2)  ! [xpi,ypi]; [xx,yy]
	! pti is nodal points in essential boundary.
	real,intent(out)::Nfe
	
	real::dist, dtemp
	integer::i
	
	!Initialization
	Nfe = dble(0)	
	
	dist = sqrt((pti(1) - pt(1))**2 + ((pti(2) - pt(2))**2))
	dtemp = max(dx,dy)	
	
	if (pti(2) - y_start < dy/10. .and. dist < dtemp) then ! pti at bottom boundary
		if (pti(1) - x_start < dx/10.) then 	! the first point
			if (pti(1) < pt(1)  .and. pt(1) < pti(1) + dx ) then ! pt falls into bottom region
				Nfe = ((pti(1) + dx) - pt(1))/(dx)
			end if
			
			if (pti(2) < pt(2) .and. pt(2) < pti(2) + dy ) then ! pt falls into left region
				Nfe = ((pti(2) + dy) - pt(2))/(dy)
			end if
			
		else if (x_end - pti(1) < dx/10.) then 	! the last point
			if (pti(1) - dx < pt(1) .and. pt(1) < pti(1)) then ! pt falls into bottom region
				Nfe = (pt(1) - (pti(1) - dx))/(dx)
			end if
			
			if (pti(2) < pt(2) .and. pt(2) < pti(2) + dy ) then ! pt falls into right region
				Nfe = ((pti(2) + dy) - pt(2))/(dy)
			end if
			
		else ! Pti in the middle 
			if (pti(1) < pt(1) .and. pt(1) < pti(1) + dx ) then  ! 
				Nfe = ((pti(1) + dx) - pt(1))/(dx)
			else
				Nfe = (pt(1) - (pti(1) - dx))/(dx)
			end if
		end if
		
	else if (y_end - pti(2) < dy/10. .and. dist < dtemp) then ! up boundary
		if (pti(1) - x_start < dx/10.) then 	! the first point
			if (pti(1) < pt(1) .and. pt(1) < pti(1) + dx ) then ! pt falls into up region
				Nfe = ((pti(1) + dx) - pt(1))/(dx)
			end if
			
			if (pti(2) - dy < pt(2) .and. pt(2) < pti(2)) then ! pt falls into left region
				Nfe = (pt(2) - (pti(2) - dy))/(dy)
			end if
			
		else if (x_end - pti(1) < dx/10.) then 	! the last point
			if (pti(1) - dx < pt(1) .and. pt(1) < pti(1)) then ! pt falls into up region
				Nfe = (pt(1) - (pti(1) - dx))/(dx)
			end if
			
			if (pti(2) - dy < pt(2) .and. pt(2) < pti(2)) then ! pt falls into right region
				Nfe = (pt(2) - (pti(2) - dy))/(dy)
			end if
			
		else ! Pti in the middle 
			if (pti(1) < pt(1) .and. pt(1) < pti(1) + dx ) then  ! 
				Nfe = ((pti(1) + dx) - pt(1))/(dx)
			else
				Nfe = (pt(1) - (pti(1) - dx))/(dx)
			end if
		end if
			
	else if (x_end - pti(1) < dx/10. .and. dist < dtemp) then ! right boundary and not include the first node
		if (pti(2) < pt(2) .and. pt(2) < pti(2) + dy ) then  ! ! Pti in the middle 
			Nfe = ((pti(2) + dy) - pt(2))/(dy)
		else
			Nfe = (pt(2) - (pti(2) - dy))/(dy)
		end if	
		
	else if (pti(1) - x_start < dx/10. .and. dist < dtemp) then ! left boundary
		if (pti(2) < pt(2) .and. pt(2) < pti(2) + dy ) then  ! ! Pti in the middle
			Nfe = ((pti(2) + dy) - pt(2))/(dy)
		else
			Nfe = (pt(2) - (pti(2) - dy))/(dy)
		end if
		
	end if
	
end subroutine

!*******************************************************************************