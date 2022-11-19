subroutine readSolData(xynodes,xyGaussPts,xyGaussWts,solU,ExactsolU)
!	this is to read in the parameters: xynodes, xyGaussPts, xyGaussWts, EBCnodes, NBCnodes ...
	use constants
	implicit none
	
	! parameters declaration
	real,intent(out)::xynodes(np,2),xyGaussPts(npGauss,2),xyGaussWts(npGauss),solU(2*np),ExactsolU(2*np)
	
	integer::i
	
!*******************************************************************************************	
! Read data from previous unit
!*******************************************************************************************	
	open(1,file="./data_pre/xynodes.txt",status="old")
	read(1,*) !the file name
	read(1,*) !the readme
	read(1,*) !the size
	do i=1,np 
		read(1,"(2E12.4)")xynodes(i,1),xynodes(i,2)
	end do
	close(1)
!*******************************************************************************************	
	open(2,file="./data_pre/xyGaussPts.txt",status="old")
	read(2,*) !the file name
	read(2,*) !the readme
	read(2,*) !the size
	do i=1,npGauss
		read(2,"(2E12.4)")xyGaussPts(i,:)
	end do
	close(2)
!*******************************************************************************************	
	open(3,file="./data_pre/xyGaussWts.txt",status="old")
	read(3,*) !the file name
	read(3,*) !the readme
	read(3,*) !the size
	do i=1,npGauss
		read(3,"(2E12.4)")xyGaussWts(i)
	end do
	close(3)
!*******************************************************************************************	
	open(15,file="./data_pro/solU.txt",status="old")
	read(15,*) !the file name
	read(15,*) !the readme
	read(15,*) !the size
	do i=1,2*np 
		read(15,"(E12.4)")solU(i)
	end do
	close(15)
!*******************************************************************************************	
	open(16,file="./data_pro/Vec_uBefo(exact).txt",status="old")
	read(16,*) !the file name
	read(16,*) !the readme
	read(16,*) !the size
	do i=1,2*np 
		read(16,"(E12.4)")ExactsolU(i)
	end do
	close(16)
	
end subroutine
