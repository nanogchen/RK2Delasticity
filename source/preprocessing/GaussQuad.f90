!*******************************************************************************

subroutine GaussQuad1D(n,GaussPtx,GaussWt)
!	This is to find the 1D standard Gauss points in [-1,1]
	! n point Gauss quadrature

	implicit none

	integer,intent(in)::n  ! n Gauss points 
	real,intent(out)::GaussPtx(n),GaussWt(n)
	integer::i
	
	!Initialization
	GaussPtx = dble(0)
	GaussWt = dble(0)
	
	if (n == 1) then
		GaussPtx(1) = dble(0)
		
		GaussWt(1) = 2.0
	
	else if (n == 2) then
		GaussPtx(1) = - 0.577350269  ! t1  ! +- sqrt(1/3)
		GaussPtx(2) = - GaussPtx(1)  ! t2
		
		GaussWt(1) = 1.0
		GaussWt(2) = 1.0
		
	else if (n == 3) then
		GaussPtx(1) = - 0.774596669   ! t1 ! +- sqrt(3/5)
		GaussPtx(2) = dble(0)  		  ! t2
		GaussPtx(3) = - GaussPtx(1)   ! t3
		
		GaussWt(1) = (5./9.)
		GaussWt(2) = (8./9.)
		GaussWt(3) = GaussWt(1)
		
	else if (n == 4) then
		GaussPtx(1) = - 0.861136312   ! t1
		GaussPtx(2) = - 0.339981044   ! t2
		GaussPtx(3) = - GaussPtx(2)   ! t3
		GaussPtx(4) = - GaussPtx(1)   ! t4
		
		GaussWt(1) = (0.347854845)
		GaussWt(2) = (0.652145155)
		GaussWt(3) = GaussWt(2)
		GaussWt(4) = GaussWt(1)
		
	else if (n == 5) then
		GaussPtx(1) = - 0.906179846   ! t1
		GaussPtx(2) = - 0.538469310   ! t2
		GaussPtx(3) = dble(0)  		  ! t3
		GaussPtx(4) = - GaussPtx(2)   ! t4
		GaussPtx(5) = - GaussPtx(1)   ! t5
		
		GaussWt(1) = (0.236926885)
		GaussWt(2) = (0.478628670)
		GaussWt(3) = (0.568888889)
		GaussWt(4) = GaussWt(2)
		GaussWt(5) = GaussWt(1)
			
	else if (n == 6) then
		GaussPtx(1) = - 0.932469514   ! t1
		GaussPtx(2) = - 0.661209385   ! t2
		GaussPtx(3) = - 0.238619186	  ! t3
		GaussPtx(4) = - GaussPtx(3)   ! t4
		GaussPtx(5) = - GaussPtx(2)   ! t5
		GaussPtx(6) = - GaussPtx(1)   ! t6
		
		GaussWt(1) = (0.171324492)
		GaussWt(2) = (0.360761573)
		GaussWt(3) = (0.467913935)
		GaussWt(4) = GaussWt(3)
		GaussWt(5) = GaussWt(2)
		GaussWt(6) = GaussWt(1)
		
	else if (n == 7) then
		GaussPtx(1) = - 0.949107912   ! t1
		GaussPtx(2) = - 0.741531186   ! t2
		GaussPtx(3) = - 0.405845151	  ! t3
		GaussPtx(4) = dble(0)   	  ! t4
		GaussPtx(5) = - GaussPtx(3)   ! t5
		GaussPtx(6) = - GaussPtx(2)   ! t6
		GaussPtx(7) = - GaussPtx(1)   ! t7
		
		GaussWt(1) = (0.129484966)
		GaussWt(2) = (0.279705391)
		GaussWt(3) = (0.381830051)
		GaussWt(4) = (0.417959183)
		GaussWt(5) = GaussWt(3)
		GaussWt(6) = GaussWt(2)
		GaussWt(7) = GaussWt(1)
		
	else if (n == 8) then
		GaussPtx(1) = - 0.960289856   ! t1
		GaussPtx(2) = - 0.796666477   ! t2
		GaussPtx(3) = - 0.525532410	  ! t3
		GaussPtx(4) = - 0.183434642	  ! t4
		GaussPtx(5) = - GaussPtx(4)   ! t5
		GaussPtx(6) = - GaussPtx(3)   ! t6
		GaussPtx(7) = - GaussPtx(2)   ! t7
		GaussPtx(8) = - GaussPtx(1)   ! t8
		
		GaussWt(1) = (0.101228536)
		GaussWt(2) = (0.222381034)
		GaussWt(3) = (0.313706646)
		GaussWt(4) = (0.362683783)
		GaussWt(5) = GaussWt(4)
		GaussWt(6) = GaussWt(3)
		GaussWt(7) = GaussWt(2)
		GaussWt(8) = GaussWt(1)
		
	else 
		write(*,*)"Not applicable yet! Exiting..."
		stop
		
	end if
	
end subroutine

!*******************************************************************************

subroutine GaussQuad2D(Pt1,Pt2,Pt3,Pt4,nGx,nGy,GaussPtx,GaussWt)
	! This is to find the 2D Gauss points in [-1,1]*[-1,1] and mapping to the real coodinate sys
	! nGx*nGy Gauss points
	! Pt1,Pt2,Pt3,Pt4 original points
	implicit none
	
	integer,intent(in)::nGx,nGy
	real,intent(in)::Pt1(2),Pt2(2),Pt3(2),Pt4(2)
	real,intent(out)::GaussPtx(nGx*nGy,2),GaussWt(nGx*nGy)  !coord in original system, not standard [-1,1] system
	
	real::xGaussPt(nGx),yGaussPt(nGy),xGaussWt(nGx),yGaussWt(nGy),Jacobian,Pt_temp(2)
	integer::i,j,k
	
	!Initialization
	GaussPtx = dble(0)
	GaussWt = dble(0)
	xGaussPt = dble(0)
	xGaussWt = dble(0)
	yGaussPt = dble(0)
	yGaussWt = dble(0)
	Jacobian = dble(0)
	Pt_temp = dble(0)
	
	call GaussQuad1D(nGx,xGaussPt,xGaussWt)
	call GaussQuad1D(nGy,yGaussPt,yGaussWt)
	
	! find the Gauss points in standard coord sys a-b
	do i=1,nGy
		do j=1,nGx
			k = (i-1)*nGx+j
			GaussPtx(k,1) = xGaussPt(j)
			GaussPtx(k,2) = yGaussPt(i)
			
			GaussWt(k) = xGaussWt(j)*yGaussWt(i)   ! w = wi*wj
		end do 
	end do	
	
	! mapping to the original x-y cood sys
	do i=1,nGy
		do j=1,nGx
			k = (i-1)*nGx+j
			call GaussMapping(Pt1,Pt2,Pt3,Pt4,Pt_temp,GaussPtx(k,:),Jacobian)
			GaussPtx(k,:) = Pt_temp				
			GaussWt(k) = GaussWt(k)*Jacobian ! w = wi*wj*J
		end do
	end do 
	
end subroutine

!*******************************************************************************

subroutine GaussMapping(Pt1,Pt2,Pt3,Pt4,Pt,PtSt,Jacobian)
! arbitrary x-y coord sys ---> a-b coord sys. 2D mapping
	implicit none
	
	real,intent(in)::PtSt(2),Pt1(2),Pt2(2),Pt3(2),Pt4(2)
!	Pt1-Pt4 are the arbitrary quadrilateral points in x-y coord sys, Pt is the Gauss pt in x-y coord sys 
	real,intent(out)::Pt(2),Jacobian ! PtSt is the mapping standard point in a-b coord sys == [-1,1]*[-1,1]

	real::dxa,dya,dxb,dyb   ! to calculate Jacobian 
	real::N1,N2,N3,N4  ! shape function, ************which is in the a-b coord sys************

	N1 = dble(0)
	N2 = dble(0)
	N3 = dble(0)
	N4 = dble(0)	
	Jacobian = dble(0)
	dxa = dble(0)
	dya = dble(0)
	dxb = dble(0)
	dyb = dble(0)
	Pt = (/ dble(0), dble(0) /)

	! shape function of the standard coodinate system
	N1 = 1./4.*(1. - PtSt(1))*(1. - PtSt(2))
	N2 = 1./4.*(1. + PtSt(1))*(1. - PtSt(2))
	N3 = 1./4.*(1. + PtSt(1))*(1. + PtSt(2))
	N4 = 1./4.*(1. - PtSt(1))*(1. + PtSt(2))

	Pt(1) = N1*Pt1(1) + N2*Pt2(1) + N3*Pt3(1) + N4*Pt4(1)
	Pt(2) = N1*Pt1(2) + N2*Pt2(2) + N3*Pt3(2) + N4*Pt4(2)

	dxa = Pt1(1)*(-1./4.*(1-PtSt(2))) + Pt2(1)*(1./4.*(1-PtSt(2))) + Pt3(1)*(1./4.*(1+PtSt(2))) + Pt4(1)*(-1./4.*(1+PtSt(2)))
	dya = Pt1(2)*(-1./4.*(1-PtSt(2))) + Pt2(2)*(1./4.*(1-PtSt(2))) + Pt3(2)*(1./4.*(1+PtSt(2))) + Pt4(2)*(-1./4.*(1+PtSt(2)))
	dxb = Pt1(1)*(-1./4.*(1-PtSt(1))) + Pt2(1)*(-1./4.*(1+PtSt(1))) + Pt3(1)*(1./4.*(1+PtSt(1))) + Pt4(1)*(1./4.*(1-PtSt(1)))
	dyb = Pt1(2)*(-1./4.*(1-PtSt(1))) + Pt2(2)*(-1./4.*(1+PtSt(1))) + Pt3(2)*(1./4.*(1+PtSt(1))) + Pt4(2)*(1./4.*(1-PtSt(1)))

	Jacobian = abs(dxa*dyb - dxb*dya)

end subroutine

!*******************************************************************************

subroutine generateGaussPts1D(np,xnodes,nG,xGaussPts,xGaussWts)
	implicit none
! 	find the gauss points in 1D problem.
	
	integer,intent(in)::np,nG ! nG is the num of Gauss points
	real,intent(in)::xnodes(np)
	real,intent(out)::xGaussPts(nG*(np-1)),xGaussWts(nG*(np-1))
	
	real::GaussPtx(nG),GaussWt(nG)
	integer::i,j
	
	!initialization
	xGaussPts = dble(0)
	xGaussWts = dble(0)
	GaussPtx = dble(0)
	GaussWt = dble(0)
	
	do i = 1,np-1          ! use nG point Gauss quadrature
		call GaussQuad1D(nG,GaussPtx,GaussWt)
		do j=1,nG
			xGaussPts(nG*(i-1)+j) = GaussPtx(j)*(xnodes(i+1)-xnodes(i))/2.0 + (xnodes(i+1)+xnodes(i))/2.
			xGaussWts(nG*(i-1)+j) = GaussWt(j)*abs((xnodes(i+1)-xnodes(i))/2.0) ! mapping back
		end do
	end do
	
end subroutine

!*******************************************************************************

subroutine generateGaussPts2D(npx,npy,xynodes,nGx,nGy,xyGaussPts,xyGaussWts)
	implicit none

	integer,intent(in)::npx,npy,nGx,nGy
	real,intent(in)::xynodes(npx*npy,2)
	real,intent(out)::xyGaussPts(nGx*nGy*(npx-1)*(npy-1),2),xyGaussWts(nGx*nGy*(npx-1)*(npy-1))

	real::GaussPts(nGx*nGy,2),GaussWt(nGx*nGy)
	integer::i,j,k,l,i_local,j_local,k_local
	integer::npxGauss ! num of Gauss points in x-axis
	
	!initialization
	xyGaussPts = dble(0)
	xyGaussWts = dble(0)
	GaussPts = dble(0)
	GaussWt = dble(0)
	npxGauss = nGx*(npx-1)
	
	write(*,*)"Generating the Gauss points..."

	do i = 1,npy-1 		! row of xynodes
		do j = 1,npx-1 		! column of xynodes
			k = (i-1)*npx + j   !!!!! index point in xynodes
			call GaussQuad2D(xynodes(k,:),xynodes(k+1,:),xynodes(k+npx+1,:),xynodes(k+npx,:),nGx,nGy,GaussPts,GaussWt)
			! Notice the order of the 4 nodal points, don't flip the order
						
			l = (i-1)*nGy*npxGauss + (j-1)*nGx + 1  !!!!! index point in xyGaussPts			
			do i_local=1,nGy
				do j_local=1,nGx
					k_local = (i_local-1)*nGx + j_local
					
					xyGaussPts(l + j_local-1 + (i_local-1)*npxGauss,:) = GaussPts(k_local,:)
					xyGaussWts(l + j_local-1 + (i_local-1)*npxGauss) = GaussWt(k_local)
				end do
			end do
		end do 
	end do 

end subroutine

!*******************************************************************************