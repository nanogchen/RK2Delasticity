!*******************************************************************************

subroutine CalExactDisp(disp, pt)
	use constants
	implicit none
	
	real,intent(in):: pt(2)
	real,intent(out):: disp(2)
	
	real::x,y
	
	!initialization
	disp = dble(0)
	x = pt(1)
	y = pt(2)
	
!!!!! for the patch test
	! linear
	! disp(1) = 0.2*x + 0.3*y
	! disp(2) = 0.1*x + 0.4*y
	! quadratic
	! disp(1) = 0.12*x + 0.14*y + 0.16*x**2 + 0.18*x*y + 0.20*y**2
	! disp(2) = 0.11*x + 0.13*y + 0.15*x**2 + 0.17*x*y + 0.19*y**2	
	
!!!! for cantilever beam: Timoshenko
	disp(1) = Pforce*y/(6.0*Emod*Imoment)*((6.0*Length-3*x)*x + (2+nu)*(y**2 - Height**2/4))
	disp(2) = -Pforce/(6.0*Emod*Imoment)*(3*nu*y**2*(Length-x) + (4+5*nu)*Height**2*x/4 + (3*Length - x)*x**2)
	
!!!! for slender beam: Euler
	! fixed EBC
	! disp(1) = dble(0)
	! disp(2) = Pforce*x**2/(6*Emod*Imoment)*(3*Length -x)

	
!!!!! for sin EBC
	! disp(1) = sin(pi*x)*sin(pi*y)
	! disp(2) = sin(pi*x)*sin(pi*y)
	
end subroutine

!*******************************************************************************

subroutine CalExactDomainDisp(npts,disp, points)
	use constants
	implicit none
	
	integer,intent(in)::npts
	real,intent(in):: points(npts,2)
	real,intent(out):: disp(2*npts)
	
	real::disp_temp(2)
	integer::i
	
	!initialization
	disp = dble(0)
	disp_temp = dble(0)
	
	do i = 1,npts
		call CalExactDisp(disp_temp, points(i,:))
		disp(2*(i-1)+1) = disp_temp(1)
		disp(2*(i-1)+2) = disp_temp(2)
	end do 
	
	! for cantilever beam	
	
end subroutine

!*******************************************************************************

subroutine CalBodyForce(force,pt)
	use constants
	implicit none
	
	real,intent(in):: pt(2)
	real,intent(out):: force(2)
	
	!initialization
	force = dble(0)
	
!!!!! quadratic
	!plane stress
	! force(1) = -51730.77
	! force(2) = -56923.08
	!plane strain
	! force(1) = -62692
	! force(2) = -66154
	
	! describe the form of the body force
	
end subroutine


!*******************************************************************************

subroutine CalTraction(force,Pt)
	use constants
!	THIS IS TO calculate the traction force on the natural boundary
	implicit none
	
	real,intent(in)::Pt(2) ! Nodal point coordinate 
	real,intent(out)::force(2)
	
	!Initialization
	force = dble(0)	

	!force = some other function (cantilever beam: Timoshenko)
	force(1) = dble(0)
	force(2) = -Pforce/(2.0*Imoment)*(Height**2/4. - Pt(2)**2)
	
	!force = cantilever beam: Euler slender beam
	! force(1) = dble(0)
	! force(2) = Pforce; !Area=1
	! force(1) = dble(0)
	! force(2) = dble(0) !use tip load

	
end subroutine

!*******************************************************************************

subroutine AddTipLoad(index, load, Vec_f)
!	this is to add the tip load boundary condition
	use constants,only:np
	implicit none

	integer,intent(in)::index
	real,intent(in)::load(2)
	real,intent(inout)::Vec_f(2*np)

	Vec_f(2*index-1:2*index) = Vec_f(2*index-1:2*index) + load

end


!*******************************************************************************

subroutine CalExactStress(Pt, stress)
	!to calulate the exact stress at a given point
	use constants
	implicit none

	real,intent(in)::Pt(2)
	real,intent(out)::stress(3)

	real::x,y

	!initialization
	stress = dble(0)
	x = pt(1)
	y = pt(2)

	stress(1) = Pforce*(Length - x)*y/Imoment
	stress(2) = dble(0);
	stress(3) = -Pforce/(2.0*Imoment)*(Height**2/4.0 - y**2)

end subroutine

!*******************************************************************************

subroutine CalExactDomainStress(npGauss, xyGaussPts, StressGauss)
	implicit none

	integer,intent(in)::npGauss
	real,intent(in)::xyGaussPts(npGauss,2)
	real,intent(out)::StressGauss(3,npGauss)

	integer::i
	real::stress_temp(3)

	!initialization
	StressGauss = dble(0)
	stress_temp = dble(0)

	do i = 1,npGauss
		call CalExactStress(xyGaussPts(i,:), stress_temp)
		StressGauss(:,i) = stress_temp
	end do


end subroutine