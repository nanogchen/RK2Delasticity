subroutine readInputData(xynodes,xyGaussPts,xyGaussWts,EBCnodes,NBCnodes,&
						EBCGaussPts,NBCGaussPts,EBCGaussWts,NBCGaussWts,OrderedNodeList)
!	this is to read in the parameters: xynodes, xyGaussPts, xyGaussWts, EBCnodes, NBCnodes ...
	use constants
	implicit none
	
	! parameters declaration
	real,intent(out)::xynodes(np,2),xyGaussPts(npGauss,2),xyGaussWts(npGauss)
	real,intent(out)::EBCnodes(npEBC,2),NBCnodes(npNBC,2)
	real,intent(out)::EBCGaussPts(npGaussEBC,2),NBCGaussPts(npGaussNBC,2)
	real,intent(out)::EBCGaussWts(npGaussEBC),NBCGaussWts(npGaussNBC)
	integer,intent(out)::OrderedNodeList(np)
	
	integer::i
	
!*******************************************************************************************	
! Read geo_data from proprocessed
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
	open(4,file="./data_pre/EBCnodes.txt",status="old")
	read(4,*) !the file name
	read(4,*) !the readme
	read(4,*) !the size
	do i=1,npEBC
		read(4,"(2E12.4)")EBCnodes(i,:)
	end do
	close(4)
!*******************************************************************************************
	open(5,file="./data_pre/NBCnodes.txt",status="old")
	read(5,*) !the file name
	read(5,*) !the readme
	read(5,*) !the size
	do i=1,npNBC
		read(5,"(2E12.4)")NBCnodes(i,:)
	end do
	close(5)
!*******************************************************************************************
	open(66,file="./data_pre/EBCGaussPts.txt",status="old")
	read(66,*) !the file name
	read(66,*) !the readme
	read(66,*) !the size
	do i=1,npGaussEBC
		read(66,"(2E12.4)")EBCGaussPts(i,:)
	end do
	close(66)
!*******************************************************************************************
	open(7,file="./data_pre/EBCGaussWts.txt",status="old")
	read(7,*) !the file name
	read(7,*) !the readme
	read(7,*) !the size
	do i=1,npGaussEBC
		read(7,"(2E12.4)")EBCGaussWts(i)
	end do
	close(7)
!*******************************************************************************************
	open(8,file="./data_pre/NBCGaussPts.txt",status="old")
	read(8,*) !the file name
	read(8,*) !the readme
	read(8,*) !the size
	do i=1,npGaussNBC
		read(8,"(2E12.4)")NBCGaussPts(i,:)
	end do
	close(8)
!*******************************************************************************************
	open(9,file="./data_pre/NBCGaussWts.txt",status="old")
	read(9,*) !the file name
	read(9,*) !the readme
	read(9,*) !the size
	do i=1,npGaussNBC
		read(9,"(2E12.4)")NBCGaussWts(i)
	end do	
	close(9)
	
!*******************************************************************************************
	open(123,file="./data_pre/OrderedNodeList.txt",status="old")
	read(123,*) !the file name
	read(123,*) !the readme
	read(123,*) !the size
	do i=1,np
		read(123,*)OrderedNodeList(i)
	end do	
	close(123)
	
	
end subroutine
