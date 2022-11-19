subroutine rearrangeMat(np,NodeList,matBefo,FBefo,Vec_UBefo,matAfter,FAfter,Vec_UAfter,Trans_matRight)
!	this is to rearrange the stiffness matrix to adjust it to block matrix
!	matBefo,FBefo,Vec_UBefo: is the K,f and u respectively, so that Ku=f and rearrange it to after-form
	implicit none 
	
	integer,intent(in)::np,NodeList(np)
	real,intent(in)::matBefo(2*np,2*np),FBefo(2*np),Vec_UBefo(2*np)
	real,intent(out)::matAfter(2*np,2*np),FAfter(2*np),Vec_UAfter(2*np)
	integer,intent(out)::Trans_matRight(2*np,2*np)
	
	integer::Trans_matLeft(2*np,2*np)  ! unit mat to arrange the origin matrix
	integer::i,j
	
	!Initialization
	Trans_matLeft = 0  ! row transformation
	Trans_matRight = 0 ! colomn transformation
	matAfter = dble(0)
	FAfter = dble(0)
	Vec_UAfter = dble(0)
	
	write(*,*)"Rearranging the system matrix and vector to impose the BC..."
	
	! To get the unit matrix to move specific row and column of the matrix 
	do i=1,np 
		do j=1,np 
			if (j == NodeList(i)) then
				Trans_matLeft(2*(i-1)+1,2*(j-1)+1) = 1
				Trans_matLeft(2*(i-1)+2,2*(j-1)+2) = 1
			end if 
		end do 
	end do
	! call saveIntData2D(1111, "Trans_matLeft.txt", 2*np, 2*np, Trans_matLeft)
	
	Trans_matRight = transpose(Trans_matLeft)
	matAfter = matmul(Trans_matLeft,matBefo)
	matAfter = matmul(matAfter,Trans_matRight)
	
	FAfter = matmul(Trans_matLeft,FBefo)
	Vec_UAfter = matmul(Trans_matLeft,Vec_UBefo)
	
	write(*,*)"Done!"
	
end subroutine