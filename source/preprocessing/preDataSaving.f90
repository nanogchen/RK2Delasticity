!*******************************************************************************

subroutine saveRealData1D(num, nameoftxt, n, xpts)
!	this is for saving the data
    implicit none

    integer,intent(in)::num, n! for file unit
    character(*),intent(in):: nameoftxt
    real,intent(in)::xpts(n)  !x-axis coordinate; , ypts(i)
	character(len=100):: Filename
	integer::i  

    write(*,*)"Saving data... ",nameoftxt

    Filename = "./data_pre/"//nameoftxt
    open(num, file=Filename)
	write(num,*)"#",nameoftxt
	write(num,*)"#SizeofTheFile"
	write(num,"(I8)")n
    do i=1, n
        write(num,"(E12.4)")xpts(i)
!         write(num,"(E12.4)",advance="no")xpts(i)
    end do
    close(num)

    write(*,*)"Saved!"

end subroutine

!*******************************************************************************

subroutine saveCoord2D(num, nameoftxt, nrow, mat)
!	this is for saving the nrow * ncol data (matrix)
    implicit none

    integer,intent(in)::num, nrow ! for file unit
    character(*),intent(in):: nameoftxt
	real,intent(in)::mat(nrow,2)  !data of matrix mat
    character(len=100):: Filename
	integer::i, j 

    write(*,*)"Saving data... ",nameoftxt

	Filename = "./data_pre/"//nameoftxt
    open(num, file=Filename)
	write(num,*)"#",nameoftxt
	write(num,*)"#SizeofTheFile"
	write(num,"(I8,I8)")nrow,2
    do i = 1, nrow
        write(num,"(2E12.4)")mat(i,:)
    end do

    close(num)
    write(*,*)"Saved!"

end subroutine

!*******************************************************************************

subroutine vtkWriteNodes(np,xynodes) ! Nodal coordinate
	implicit none
	
	integer,intent(in)::np
	real,intent(in)::xynodes(np,2)
	
	character(len=100):: Filename
	integer::num,i
	
	Filename = "./data_pre/Geometry.vtk"
	num = 111
	
	open(num, file=Filename)
	write(num,'(A26)')"# vtk DataFile Version 4.0"
	write(num,*)"MyRK 2D file"
	write(num,'(A5)')"ASCII"
	write(num,'(A25)') 'DATASET UNSTRUCTURED_GRID'
	write(num,'(A7,I10,A10)')"POINTS ",np,"double"
	! coordinate
	do i=1, np
		write(num,"(3E12.4)")xynodes(i,1), xynodes(i,2), dble(0)
    end do
	
	! cell
	write(num,*)' '
	write(num,*)"CELLS", np, np*2
	do i=1, np
		write(num,*)1, i-1
    end do
	
	! cell types
	write(num,*)' '
	write(num,*)"CELL_TYPES", np
	do i=1, np
		write(num,*)1
    end do

    close(num)
	

end subroutine

!*******************************************************************************

subroutine vtkWritePoints(npGauss,xyGaussPts,xyGaussWts) ! Nodal coordinate
	implicit none
	
	integer,intent(in)::npGauss
	real,intent(in)::xyGaussPts(npGauss,2),xyGaussWts(npGauss)
	
	character(len=100):: Filename
	integer::num,i
	
	Filename = "./data_pre/Gauss.vtk"
	num = 222
	
	open(num, file=Filename)
	write(num,'(A26)')"# vtk DataFile Version 4.0"
	write(num,'(A12)')"MyRK 2D file"
	write(num,'(A5)')"ASCII"
	write(num,'(A25)') 'DATASET UNSTRUCTURED_GRID'
	write(num,'(A7,I10,A10)')"POINTS ",npGauss,"double"
	! coordinate
	do i=1, npGauss
		write(num,"(3E12.4)")xyGaussPts(i,1), xyGaussPts(i,2), dble(0)
    end do
	
	! cell
	write(num,*)' '
	write(num,*)"CELLS", npGauss, npGauss*2
	do i=1, npGauss
		write(num,*)1, i-1
    end do
	
	! cell types
	write(num,*)' '
	write(num,*)"CELL_TYPES", npGauss
	do i=1, npGauss
		write(num,*)1
    end do
	
	! point data
	write(num,*)' '
	WRITE(num,'(A11, I10)') 'POINT_DATA ', npGauss
	
	! Gauss weight
	write(num,*)' '
	WRITE(num,'(A8, A30, A10)')    'SCALARS ',   'gauss_weight',     'double'
	WRITE(num,'(A20)') 'LOOKUP_TABLE default'
	do i=1, npGauss
		write(num,*)xyGaussWts(i)
    end do	

    close(num)	

end subroutine


!*******************************************************************************

subroutine saveRealData2D(num, nameoftxt, nrow, ncol, mat)
!	this is for saving the nrow * ncol data (matrix)
    implicit none

    integer,intent(in)::num, nrow, ncol ! for file unit
    character(*),intent(in):: nameoftxt
	real,intent(in)::mat(nrow,ncol)  !data of matrix mat
    character(len=100):: Filename
	integer::i, j 

    write(*,*)"Saving data... ",nameoftxt

    Filename = "./data_pre/"//nameoftxt
    open(num, file=Filename)
	write(num,*)"#",nameoftxt
	write(num,*)"#SizeofTheFile"
	write(num,"(I8,I8)")nrow,ncol
    do i = 1, nrow
    	do j = 1, ncol
        	write(num,"(E12.4)", advance = "no")mat(i,j)
    	end do
    	write(num,"(/)",advance = "no")
    end do

    close(num)
    write(*,*)"Saved!"

end subroutine

!*******************************************************************************

subroutine saveIntData1D(num, nameoftxt, n, xpts)
!	this is for saving the data
    implicit none

    integer,intent(in)::num, n, xpts(n)! for file unit
    character(*),intent(in):: nameoftxt
    character(len=100):: Filename
	integer::i  

    write(*,*)"Saving data... ",nameoftxt

    Filename = "./data_pre/"//nameoftxt
    open(num, file=Filename)
	write(num,*)"#",nameoftxt
	write(num,*)"#SizeofTheFile"
	write(num,"(I8)")n
    do i=1, n
        write(num,*)xpts(i)
    end do
    close(num)

    write(*,*)"Saved!"

end subroutine

!*******************************************************************************

subroutine saveIntData2D(num, nameoftxt, nrow, ncol, mat)
!	this is for saving the nrow * ncol data (matrix)
    implicit none

    integer,intent(in)::num, nrow, ncol,mat(nrow,ncol) !data of matrix mat
    character(*),intent(in):: nameoftxt
    character(len=100):: Filename	
	integer::i, j 

    write(*,*)"Saving data... ",nameoftxt

    Filename = "./data_pre/"//nameoftxt
    open(num, file=Filename)
	write(num,*)"#",nameoftxt
	write(num,*)"#SizeofTheFile"
	write(num,"(I8,I8)")nrow,ncol
    do i = 1, nrow
    	do j = 1, ncol
        	write(num,"(I5)", advance = "no")mat(i,j)
    	end do
    	write(num,"(/)",advance = "no")
    end do

    close(num)
    write(*,*)"Saved!"

end subroutine

!*******************************************************************************
