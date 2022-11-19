program RK2Delasticity_preprocessing
	implicit none
	
interface
	subroutine generateBCNodes(np,npx,npy,xynodes,BoundaryList,npEBC,EBCnodes,npNBC,NBCnodes)
	integer,intent(in)::np,npx,npy
	real,intent(in)::xynodes(np,2)
	integer,intent(in)::BoundaryList(2*(npx+npy)-4)
	integer,intent(out)::npEBC,npNBC
    real,allocatable,intent(out)::EBCnodes(:,:),NBCnodes(:,:)
	end subroutine	
	!subroutine 2 ! for Gauss Points
end interface
	
! !*******************************************************************************
! ! all the calculations, like FindKf and L2Norm, are done at Gauss Points
! !*******************************************************************************
	real::time_begin, time_end
	real,allocatable::xynodes(:,:),xyGaussPts(:,:),xyGaussWts(:)
	integer,allocatable::NeighborListInner(:,:),NeighborListEBC(:,:),NeighborListNBC(:,:)
	integer,allocatable::BoundaryList(:)
	real,allocatable::BoundaryNodes(:,:),BoundaryGaussPts(:,:),BoundaryGaussWts(:)	
	real,allocatable::EBCnodes(:,:),NBCnodes(:,:)

!###### problem domain ######		
	real::x_start, x_end
	real::y_start, y_end
	integer::npx,npy,np 						! node numbers
	integer::nGx,nGy,npxGauss,npyGauss,npGauss 	! nGx*nGy Gauss points in a cell
	integer::nGBC,npBC,npInner 					! n-points Gauss quadrature points in the Essential boundary	
	integer::npEBC,npNBC
	integer::npGaussEBC,npGaussNBC
	integer::NumNeighborInner,NumNeighborEBC,NumNeighborNBC
	real::dilation,Emod,nu,beta,Pforce
	integer::i_model,i_algorithm
	
	call CPU_TIME(time_begin)
	
! ******************************************************************
! READ IN THE CONTROL DATA OF THE GEOMETRY: control.dat
! ******************************************************************	
	! open(222,file="control.dat",status="old")    
    ! read(222,*)x_start,x_end
    ! read(222,*)y_start,y_end
	! read(222,*)npx,npy
	! read(222,*)nGx,nGy
	! read(222,*)nGBC
	! read(222,*)dilation
	! read(222,"(E8.2)")Emod
	! read(222,*)nu
	! read(222,*)i_model
	! read(222,*)i_algorithm
	! read(222,"(E8.2)")beta

    ! close(222)
	
	! call the subroutine to read in the control parameters
	call read_control(x_start,x_end,y_start,y_end,npx,npy,nGx,nGy,nGBC,dilation,Emod,nu,Pforce,i_model,i_algorithm,beta)
	
! ******************************************************************
! to determine the other parameter
! ******************************************************************
	np = npx*npy
	
	npxGauss = nGx*(npx-1)
	npyGauss = nGy*(npy-1)
	npGauss = npxGauss*npyGauss
	
	npBC = 2*(npx+npy)-4
	npInner = np - npBC
	
	! NumNeighborInner = (( int(dilation) + 1 )*2)**2*nGx*nGy
	! NumNeighborEBC = 8*nGBC ! the maximum num of neighbors for an essential boundary node
	! NumNeighborNBC = 8*nGBC ! the maximum num of neighbors for a natural boundary node		
	
! ******************************************************************
! allocate the array
! ******************************************************************
	allocate(xynodes(np,2))
	allocate(xyGaussPts(npGauss,2))
	allocate(xyGaussWts(npGauss))
	allocate(BoundaryList(2*(npx+npy)-4))
	allocate(BoundaryNodes(2*(npx+npy)-4,2))
	allocate(BoundaryGaussPts(2*nGBC*(npx+npy-2),2))
	allocate(BoundaryGaussWts(2*nGBC*(npx+npy-2)))
	! allocate(NeighborListInner(np,NumNeighborInner))
	! allocate(NeighborListEBC(np,NumNeighborEBC))
	! allocate(NeighborListNBC(np,NumNeighborNBC))
	
! preprocessing*****************************************************************
! preprocessing*****************************************************************
! preprocessing*****************************************************************
	call generate2dPts(x_start,x_end,npx,y_start,y_end,npy,xynodes)
	call generateGaussPts2D(npx,npy,xynodes,nGx,nGy,xyGaussPts,xyGaussWts) !for inner area integral 
	call GetWholeBoundary(np,npx,npy,npInner,x_start,x_end,y_start,y_end,xynodes,BoundaryList,BoundaryNodes)
	call generateBCNodes(np,npx,npy,xynodes,BoundaryList,npEBC,EBCnodes,npNBC,NBCnodes) ! in the interface
	call generateBCGaussPts(npx,npy,nGBC,BoundaryNodes,BoundaryGaussPts,BoundaryGaussWts)	
	! call SearchNeighborList(np,npGauss,npEBC,NumNeighborInner,NumNeighborEBCxynodes,xyGaussPts,EBCGaussPts,NeighborListInner,NeighborListEBC)
	
	call saveCoord2D(1,"xynodes.txt",np,xynodes)
	call saveCoord2D(2,"xyGaussPts.txt",npGauss,xyGaussPts)
	call saveRealData1D(3,"xyGaussWts.txt",npGauss,xyGaussWts)
	call saveCoord2D(4,"EBCnodes.txt",npEBC,EBCnodes)
	call saveCoord2D(5,"NBCnodes.txt",npNBC,NBCnodes)
	
	! check whether the EBC or NBC is a closed boundary
	if (npEBC == 0) then
		! all NBC
		npGaussEBC = 0
		if (npNBC == npBC) then
			npGaussNBC = npNBC*nGBC
		else if(npNBC /= 0) then
			npGaussNBC = (npNBC - 1)*nGBC
		end if
		
	else if(npNBC == 0 .and. npEBC == npBC) then
		! all EBC
		npGaussEBC = npEBC*nGBC
	else
		npGaussEBC = (npEBC - 1)*nGBC
	end	if
	
	if (npNBC == 0) then 
		! all EBC
		npGaussNBC = 0
		if (npEBC == npBC) then
			npGaussEBC = npEBC*nGBC
		else if(npEBC /= 0) then
			npGaussEBC = (npEBC - 1)*nGBC
		end if
		
	elseif(npEBC == 0 .and. npNBC == npBC) then
		! all NBC
		npGaussNBC = npNBC*nGBC
	else
		npGaussNBC = (npNBC - 1)*nGBC
	end if

	
	! write out the parameters to be used by processing unit
	call constantFileWrite(x_start,x_end,y_start,y_end,npx,npy,nGx,nGy,nGBC,npEBC,npNBC,npGaussEBC,npGaussNBC, &
							dilation,Emod,nu,Pforce,i_model,i_algorithm,beta)
	
	call CPU_TIME(time_end)
	write(*,*) "Computation time is:", time_end-time_begin
	
	deallocate(xynodes)
	deallocate(EBCnodes)
	deallocate(NBCnodes)
	deallocate(xyGaussPts)
	deallocate(xyGaussWts)
	deallocate(BoundaryList)
	deallocate(BoundaryNodes)
	deallocate(BoundaryGaussPts)
	deallocate(BoundaryGaussWts)
	! deallocate(NeighborListInner)
	! deallocate(NeighborListEBC)
	! deallocate(NeighborListNBC)

end program

!********************************************************************************************************************

subroutine generateBCNodes(np,npx,npy,xynodes,BoundaryList,npEBC,EBCnodes,npNBC,NBCnodes)
!   this is to find the EBC NBC nodes
!   classify the nodes: essential/natural boundary nodes and inner domain nodes
    implicit none    
    
	integer,intent(in)::np,npx,npy
    real,intent(in)::xynodes(np,2)
	integer,intent(in)::BoundaryList(2*(npx+npy)-4)
	integer,intent(out)::npEBC,npNBC
    real,allocatable,intent(out)::EBCnodes(:,:),NBCnodes(:,:)
	
	integer::OrderedNodeList(np)
  	integer,allocatable::EBCList(:),NBCList(:)
    integer::array_EBC(4),array_NBC(4)       ! specify the boundary condition
    integer::i,AllocateStatus
	
	! Initialization
	npEBC = 0
	npNBC = 0
	OrderedNodeList = 0
	
	write(*,*)"Specifying the EBC/NBC nodes..."   
    
! ******************************************************************
! READ IN THE BOUNDARY CONDITION: BCReadIn.dat
! ******************************************************************
    open(111,file="BCReadIn.dat",status="old")    
    read(111,*)array_EBC
    read(111,*)array_NBC

    do i = 1,4
        if (array_EBC(i) == 1 .and. array_NBC(i)==1) then
            stop "Fatal error, Boundary override !!!!"
        end if
    end do

    close(111)
	
! ******************************************************************
! STEP1: FIND THE BOUNDARY NODE LIST
! ******************************************************************

    call specifyBCList(npx,npy,BoundaryList,array_EBC,npEBC,EBCList) ! for EBC
    call specifyBCList(npx,npy,BoundaryList,array_NBC,npNBC,NBCList) ! for NBC

! ******************************************************************
! STEP2: FIND THE BOUNDARY NODE COORDINATES FROM THE NODE LIST
! ******************************************************************
	
    allocate( EBCnodes(1:npEBC,1:2) )	
    allocate( NBCnodes(1:npNBC,1:2) )
	
    EBCnodes(1:npEBC,:) = xynodes(EBCList(:),:)
    NBCnodes(1:npNBC,:) = xynodes(NBCList(:),:)	
	
	write(*,*)"Done of specifyBCNodes!" 
	
!************************************************************
!*********** save for transformation method usage ***********
!************************************************************
	call getOrderedNodeList(npx,npy,npEBC,EBCList,OrderedNodeList)
	call saveIntData1D(123, "OrderedNodeList.txt", np, OrderedNodeList)

! deallocate the dynamic array
    deallocate(EBCList)
    deallocate(NBCList)

    ! deallocate(EBCnodes)
    ! deallocate(NBCnodes)
!THIS DEALLOCATE CANNOT BE USED SINCE IT IS PASSED OUT ARRAY

contains

	subroutine specifyBCList(npx,npy,BoundaryList,array_BC,npBC_local,BCList)
!   this is to find the EBC or NBC number list
		implicit none
		
		integer,intent(in)::npx,npy
		integer,intent(in)::BoundaryList(2*(npx+npy)-4),array_BC(4)
		integer,intent(out)::npBC_local
		integer,allocatable,intent(out)::BCList(:)		

		! start and end coordinate of the 4 boundaries
		integer::indexVector(4,2)			 
		integer::tempList(2*(npx+npy)-4),tempList2(2*(npx+npy)-4)
		integer::i,counter
		
		! Initialization
		! Before this code, BoundaryList should obtained from other subroutine
		npBC_local = 0
		counter = 0
		tempList = BoundaryList
		indexVector=reshape((/1, npx, npx+npy-1, 2*npx+npy-2, npx, npx+npy-1, 2*npx+npy-2, 2*(npx+npy)-4/), (/4,2/))

!*************************************************
! step1: scheme A
!*************************************************
		do i=1,4
			if (array_BC(i) == 0) then
				tempList(indexVector(i,1)+1:indexVector(i,2)-1) = 0 ! the middle elements are set to be 0

				if (i==1) then
					tempList(indexVector(i,1)) = 0                              ! set the left corner point
					if (array_BC(i+1) == 0) tempList(indexVector(i,2)) = 0     ! set the right corner point
					cycle
				end if
				
				if (i==4) then
					tempList(indexVector(i,2)) = 0                              ! set the bottom corner point
					if (array_BC(i-1) == 0) tempList(indexVector(i,1)) = 0     ! set the up corner point
					cycle
				end if

				if (array_BC(i-1) == 0) tempList(indexVector(i,1)) = 0     ! set the left corner point
				if (array_BC(i+1) == 0) tempList(indexVector(i,2)) = 0     ! set the left corner point
				
			end if 
		end do

		! to get the length of the BCList
		tempList2 = tempList  
		where (tempList2 /= 0)
			tempList2 = 1
		end where
		
		npBC_local = sum(tempList2)
		
		if (array_BC(4) == 1 .and. array_BC(1) == 0) then
			npBC_local = npBC_local + 1
		end if
		
		allocate(BCList(npBC_local))

	!     write(*,*)tempList

!*************************************************
! step2: scheme B
!*************************************************
		do i=1,2*(npx+npy)-4 
			if (tempList(i) /= 0) then
				counter = counter + 1
				BCList(counter) = tempList(i)
			end if 
		end do
		
		if (array_BC(4) == 1 .and. array_BC(1) == 0) then
			counter = counter + 1
			BCList(counter) = BoundaryList(1)
		end if

		! judge
		if (counter /= npBC_local) then
			stop "Something wrong happens in the boundary category..."
		end if

		! write(*,*)BCList
		
	end subroutine specifyBCList

end subroutine
