!*******************************************************************************

subroutine CalApproxValue_transformation(np,xynodes,InvMatrix_N,solU,npAppro,xypoints,uapprox)
!	this is to calculate the approx value by transformation method
	implicit none
	
	integer,intent(in)::np,npAppro
	real,intent(in)::InvMatrix_N(2*np,2*np),xynodes(np,2),solU(2*np), xypoints(npAppro,2)
	real,intent(out)::uapprox(2*npAppro)   ! approx points
	
	real::N,Ndxy(2),Ni(np),D(2*np),phi_mat(2,2)
	integer::i,j
	
	write(*,*)"Calculating approx value..."
	
	!initialization
	! Ni = dble(0)	! this is the [N1(x),N2(x),...,Nn(x)] vector
	D = dble(0)	! this is the nodal value, u(x)=sum(Ni*D)
	N = dble(0)
	Ndxy = dble(0)
	uapprox = dble(0) ! approx value at xypoints(i,:)
	
	D = matmul(InvMatrix_N,solU)
	
	do j=1,npAppro   ! node value
		do i = 1,np 
			call CalN(N,Ndxy,xynodes,xynodes(i,1),xynodes(i,2),xypoints(j,1),xypoints(j,2))
			phi_mat = reshape((/N,0.,0.,N/),(/2,2/))
			
			uapprox(2*(j-1)+1:2*(j-1)+2) = uapprox(2*(j-1)+1:2*(j-1)+2) + matmul(phi_mat,D(2*(i-1)+1:2*(i-1)+2))
		end do
	end do	
	
	write(*,*)"Done!"

end subroutine

!*******************************************************************************

subroutine CalApproxValue(np,xynodes,solU,npAppro,xypoints,uapprox)
!	this is to calculate the approx value by nodal approximation N=sum(NI*UI)
	implicit none
	
	integer,intent(in)::np,npAppro
	real,intent(in)::xynodes(np,2),solU(2*np), xypoints(npAppro,2)
	real,intent(out)::uapprox(2*npAppro)   ! approx points
	
	real::N,Ndxy(2),phi_mat(2,2)
	integer::i,j
	
	!initialization
	N = dble(0)
	Ndxy = dble(0)
	phi_mat = dble(0)
	uapprox = dble(0)
	
	do j=1,npAppro   ! node value
		do i = 1,np
			call CalN(N,Ndxy,xynodes,xynodes(i,1),xynodes(i,2),xypoints(j,1),xypoints(j,2))
			
			phi_mat = reshape((/N,0.,0.,N/),(/2,2/))
			
			uapprox(2*(j-1)+1:2*(j-1)+2) = uapprox(2*(j-1)+1:2*(j-1)+2) + matmul(phi_mat,solU(2*(i-1)+1:2*(i-1)+2))
			
		end do
	end do	
		
end subroutine

!*******************************************************************************