subroutine generateBCGaussPts(npx,npy,nGBC,BoundaryNodes,BoundaryGaussPts,BoundaryGaussWts)
	! This is to generate Gauss points along the essential/natural boundary (OR line)
	implicit none
	
	integer,intent(in)::npx,npy,nGBC
	real,intent(in)::BoundaryNodes(2*(npx+npy)-4,2)	
	real,intent(out)::BoundaryGaussPts(2*nGBC*(npx+npy-2),2),BoundaryGaussWts(2*nGBC*(npx+npy-2))
	
	real::GaussPtx(nGBC),GaussWt(nGBC),lambda(nGBC) ! lambda is the length factor vector of each Gauss point
	integer::array_EBC(4),array_NBC(4)
	integer::i,j,k
	
	! Initialization
	lambda = dble(0)
	
	call GaussQuad1D(nGBC,GaussPtx,GaussWt)
	! change into length factor
	do i=1,nGBC
		lambda(i) = (GaussPtx(i) - (-dble(1)))/dble(2)  
	end do
	
	! for the front segments
	do i=1,(2*(npx+npy)-4)-1
		do j=1,nGBC
			k = (i-1)*nGBC + j
			BoundaryGaussPts(k,1) = (1. - lambda(j))*BoundaryNodes(i,1) + lambda(j)*BoundaryNodes(i+1,1)
			BoundaryGaussPts(k,2) = (1. - lambda(j))*BoundaryNodes(i,2) + lambda(j)*BoundaryNodes(i+1,2)
			BoundaryGaussWts = GaussWt(j)* &
			sqrt((BoundaryNodes(i+1,1)-BoundaryNodes(i,1))**2 + (BoundaryNodes(i+1,2)-BoundaryNodes(i,2))**2)/2.0
		end do 
	end do
	
	! for the last segment
	i = 2*(npx+npy)-4
	do j=1,nGBC
		k = (i-1)*nGBC + j
		BoundaryGaussPts(k,1) = (1. - lambda(j))*BoundaryNodes(i,1) + lambda(j)*BoundaryNodes(1,1)
		BoundaryGaussPts(k,2) = (1. - lambda(j))*BoundaryNodes(i,2) + lambda(j)*BoundaryNodes(1,2)
		BoundaryGaussWts = GaussWt(j)* &
			sqrt((BoundaryNodes(1,1)-BoundaryNodes(i,1))**2 + (BoundaryNodes(1,2)-BoundaryNodes(i,2))**2)/2.0
	end do
	! call saveCoord2D(3333,"BoundaryGaussPts.txt",2*nGBC*(npx+npy-2),BoundaryGaussPts)
	
! ******************************************************************
! READ IN THE BOUNDARY CONDITION: BCReadIn.dat
! ******************************************************************
    open(101,file="BCReadIn.dat",status="old")    
    read(101,*)array_EBC
    read(101,*)array_NBC
	
	do i = 1,4
        if (array_EBC(i) == 1 .and. array_NBC(i) == 1) then
            stop "Fatal error, Boundary override !!!!"
        end if
    end do

    close(111)
	
! ******************************************************************
! GENERATING ESSENTIAL/NATURAL BOUNDARY GAUSS POINTS 
! ******************************************************************	
	call specifyBCGaussPts(npx,npy,nGBC,666,7,"EBC",array_EBC,BoundaryGaussPts,BoundaryGaussWts)
	call specifyBCGaussPts(npx,npy,nGBC,8,9,"NBC",array_NBC,BoundaryGaussPts,BoundaryGaussWts)
	
end subroutine

!******************************************************************************************************************

subroutine specifyBCGaussPts(npx,npy,nGBC,num1,num2,nameoftxt,array_BC,BoundaryGaussPts,BoundaryGaussWts)
	! This is to generate Gauss points along the essential boundary (OR line)
	implicit none
	
	integer,intent(in)::npx,npy,nGBC
	integer,intent(in)::array_BC(4),num1,num2 ! use for out put
	character(*),intent(in):: nameoftxt
	real,intent(in)::BoundaryGaussPts(2*nGBC*(npx+npy-2),2),BoundaryGaussWts(2*nGBC*(npx+npy-2))
	
	character(len=100):: FilenamePts,FilenameWts
	integer::i,counter,npGBC_local ! num of EBC nodes and Gauss points quadrature used OBTAINED FROM READ-IN DATA
	integer::flag(2*nGBC*(npx+npy-2)) ! decide the particle 
	real,allocatable::BCGPts_local(:,:),BCGWts_local(:)
	integer::indexVector(4,2) 
	indexVector = reshape((/1, (npx-1)*nGBC+1, (npx-1 + npy-1)*nGBC +1, (2*(npx-1) + npy-1)*nGBC +1, &
				(npx-1)*nGBC, (npx-1 + npy-1)*nGBC, (2*(npx-1) + npy-1)*nGBC, 2*(npx-1 + npy-1)*nGBC /),(/4,2/))	
				
	write(*,*)"Begin of specifyBC Gauss points!" 
	
	! Initialization
	FilenamePts = nameoftxt//"GaussPts.txt"
	FilenameWts = nameoftxt//"GaussWts.txt"
	counter = 0
	do i=1,2*nGBC*(npx+npy-2)
		flag(i) = i ! default is included
	end do 
	
	! make a judgement about the length of the array
	npGBC_local = (array_BC(1)+array_BC(3))*(npx-1)*nGBC + (array_BC(2)+array_BC(4))*(npy-1)*nGBC
	
	allocate(BCGPts_local(npGBC_local,2))
	allocate(BCGWts_local(npGBC_local))

!*************************************************
! step1: scheme A (set non-specific GPts = 0)
!*************************************************
	do i=1,4
		if (array_BC(i) == 0) then
			flag(indexVector(i,1):indexVector(i,2)) = 0
		end if 
	end do 
	
!*************************************************
! step1: scheme B
!*************************************************	
	do i = 1,2*nGBC*(npx+npy-2)
		if (flag(i) /= 0) then 
			counter = counter + 1
			BCGPts_local(counter,:) = BoundaryGaussPts(flag(i),:)
			BCGWts_local(counter) = BoundaryGaussWts(flag(i))
		end if
	end do 
	
	if (counter /= npGBC_local) then
		stop "Something wrong happens in the EBC/NBC Gauss points generating..."
	end if
	
	call saveCoord2D(num1,FilenamePts,npGBC_local,BCGPts_local)
	call saveRealData1D(num2,FilenameWts,npGBC_local,BCGWts_local)
	
	write(*,*)"Done of specifyBC Gauss points!" 
	
	deallocate(BCGPts_local)
	deallocate(BCGWts_local)	

end subroutine

!******************************************************************************************************************