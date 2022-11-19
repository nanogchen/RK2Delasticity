!*******************************************************************************

subroutine getStrain(np, npGauss, xynodes, xyGausspts, solU, InvMatrix_N, &
					strainGauss)
	implicit none

	integer,intent(in)::np,npGauss
	real,intent(in)::InvMatrix_N(2*np,2*np),xynodes(np,2), xyGausspts(npGauss, 2), solU(2*np)
	real,intent(out)::strainGauss(3,npGauss)   ! approx points

	real::N,Ndxy(2),D(2*np),B_mat(3,2)
	integer::i,j

	!initialization
	! Ni = dble(0)	! this is the [N1(x),N2(x),...,Nn(x)] vector
	D = dble(0)	! this is the nodal value, u(x)=sum(Ni*D)
	N = dble(0)
	Ndxy = dble(0)
	B_mat = dble(0)
	strainGauss = dble(0) ! strain value at Gauss Points(i,:)
	
	D = matmul(InvMatrix_N,solU) ! to get the nodal coefficients

	do i=1,npGauss
		do j=1,np
			call CalN(N,Ndxy,xynodes,xynodes(j,1),xynodes(j,2),xyGausspts(i,1),xyGausspts(i,2))
			B_mat = reshape((/ Ndxy(1),0.0, Ndxy(2),0.0,Ndxy(2),Ndxy(1)/),(/3,2/))

			strainGauss(:,i) = strainGauss(:,i) + matmul(B_mat,D(2*(j-1)+1:2*(j-1)+2))
		end do
	end do


end 

!*******************************************************************************

subroutine getStress(npGauss, strainGauss,&
					stressGauss)
	use constants,only: i_model, HookMat_strain, HookMat_stress
	implicit none

	integer,intent(in)::npGauss
	real,intent(in)::strainGauss(3,npGauss)
	real,intent(out)::stressGauss(3,npGauss)

	integer::i
	real::HookMat(3,3)

	!initialization
	HookMat = dble(0) 
	stressGauss = dble(0)

	!to get the material tensor
	if (i_model == 1) then
		HookMat = HookMat_stress
	else if (i_model == 2) then
		HookMat = HookMat_strain
	else 
		stop "Something wrong happen in choosing the plane stress/strain model."
	end if

	do i=1,npGauss		
		stressGauss(:,i) = matmul(HookMat,strainGauss(:,i))
	end do

end