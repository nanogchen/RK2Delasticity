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

    Filename = "./data_post/"//nameoftxt
    open(num, file=Filename)
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

	Filename = "./data_post/"//nameoftxt
    open(num, file=Filename)
    do i = 1, nrow
        write(num,"(2E12.4)")mat(i,:)
    end do

    close(num)
    write(*,*)"Saved!"

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

    Filename = "./data_post/"//nameoftxt
    open(num, file=Filename)
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

    Filename = "./data_post/"//nameoftxt
    open(num, file=Filename)
    do i=1, n
        write(num,"(I6)")xpts(i)
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

    Filename = "./data_post/"//nameoftxt
    open(num, file=Filename)
    do i = 1, nrow
    	do j = 1, ncol
        	write(num,"(I6)", advance = "no")mat(i,j)
    	end do
    	write(num,"(/)",advance = "no")
    end do

    close(num)
    write(*,*)"Saved!"

end subroutine

!*******************************************************************************
