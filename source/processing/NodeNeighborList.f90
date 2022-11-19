!*******************************************************************************

subroutine SearchNeighborList(xynodes,xyGaussPts,EBCGaussPts,NBCGaussPts,NeighborListInner,NeighborListEBC,NeighborListNBC)
!	this is to find neighbor list Gauss Points for each nodes
	use constants,only:np,npGauss,npGaussEBC,npGaussNBC,NumNeighborInner,NumNeighborEBC,NumNeighborNBC,hx,hy
	implicit none 
	
	real,intent(in)::xynodes(np,2),xyGaussPts(npGauss,2),EBCGaussPts(npGaussEBC,2),NBCGaussPts(npGaussNBC,2)
	integer,intent(out)::NeighborListInner(np,NumNeighborInner),NeighborListEBC(np,NumNeighborEBC),NeighborListNBC(np,NumNeighborNBC)
	
	integer::i,j,counter
	real::distx,disty
	
	!Initialization
	NeighborListInner = 0
	NeighborListEBC = 0
	counter = 0
	distx = dble(0)
	disty = dble(0)
	
	write(*,*)"Finding neighbor nodes list..."
	
	! find neighbor list of Inner domain Gauss points for each nodes
	! for rectangular phi = phi(x)*phi(y)
	do i=1,np
		do j = 1,npGauss
			! mappint i --> j
			distx = abs(xynodes(i,1) - xyGaussPts(j,1))
			disty = abs(xynodes(i,2) - xyGaussPts(j,2))
			if (distx <= hx .and. disty <= hy) then
				!add to list NeighborListInner(i)
				counter = counter + 1
				NeighborListInner(i,counter) = j
			end if
		end do 
		
		counter = 0 ! reinitialization
	end do
	
	! find neighbor list of EBC Gauss points for each nodes
	! for rectangular phi = phi(x)*phi(y)
	do i=1,np
		do j = 1,npGaussEBC
			! mappint i --> j
			distx = abs(xynodes(i,1) - EBCGaussPts(j,1))
			disty = abs(xynodes(i,2) - EBCGaussPts(j,2))
			if (distx <= hx .and. disty <= hy) then
				!add to list NeighborListInner(i)
				counter = counter + 1
				NeighborListEBC(i,counter) = j
			end if
		end do 
		
		counter = 0 ! reinitialization
	end do
	
	! find neighbor list of NBC Gauss points for each nodes
	! for rectangular phi = phi(x)*phi(y)
	do i=1,np
		do j = 1,npGaussNBC
			! mappint i --> j
			distx = abs(xynodes(i,1) - NBCGaussPts(j,1))
			disty = abs(xynodes(i,2) - NBCGaussPts(j,2))
			if (distx <= hx .and. disty <= hy) then
				!add to list NeighborListInner(i)
				counter = counter + 1
				NeighborListNBC(i,counter) = j
			end if
		end do 
		
		counter = 0 ! reinitialization
	end do
	
	write(*,*)"Done ..."
	
end subroutine

!*******************************************************************************