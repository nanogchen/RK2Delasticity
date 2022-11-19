!*******************************************************************************

subroutine L2Norm(L2error,npGauss,xyGaussWts,uexact,uapprox)
! it is to calculate the L2-norm of 2 data: exact and approx.
	implicit none
	
	integer,intent(in)::npGauss
	real,intent(in)::xyGaussWts(npGauss),uexact(npGauss),uapprox(npGauss) ! exact and approx value at gauss points
	real,intent(out)::L2error

	integer::i,j
	
	write(*,*)"Calculating L2-norm..."
	
	L2error = dble(0)

	!	F = (uapprox - uexact)^2
	do i= 1,npGauss
		L2error = L2error + (uapprox(i)-uexact(i))**2*xyGaussWts(i)
	end do 	

	L2error = sqrt(L2error)	
	
	write(*,*)"Done!"
	
end subroutine

!*******************************************************************************