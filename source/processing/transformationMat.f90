!*******************************************************************************

subroutine transforMatrix(np,xynodes,InvMatrix_N,KK)
! this is to realize tranformation method to impose the essential boundary condition 
	implicit none
	
	integer,intent(in)::np
	real,intent(in)::xynodes(np,2)
	real,intent(out)::InvMatrix_N(2*np,2*np) !this is the transformation matrix
	real,intent(inout)::KK(2*np,2*np)  ! KK is the direct system matrix mat_K, after the subroutine, it is transformed.
	
	real::Trans_mat(2*np,2*np),N,Ndxy(2),phi_mat(2,2)
	integer::i,j 
	
	!initialization
	Trans_mat = dble(0)
	InvMatrix_N = dble(0)
	N = dble(0)
	Ndxy = dble(0)
	phi_mat = dble(0)
	
	do i=1,np
		do j=1,np 
			call CalN(N, Ndxy, xynodes, xynodes(j,1), xynodes(j,2), xynodes(i,1), xynodes(i,2))
			phi_mat = reshape((/N, 0.0, 0.0, N/),(/2,2/))
			Trans_mat(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2)= phi_mat
		end do
	end do 
	
	call inverse(Trans_mat,2*np,InvMatrix_N)
	
	KK = matmul(KK,InvMatrix_N)
	
end subroutine
