!********************************************************************************************************************

subroutine GetWholeBoundary(np,npx,npy,npInner,x_start,x_end,y_start,y_end,xynodes,BoundaryList,BoundaryNodes)
!   this is to get the boundary list of node number
    implicit none
    
	integer,intent(in)::np,npx,npy,npInner
	real,intent(in)::x_start,x_end,y_start,y_end
	real,intent(in)::xynodes(np,2)
    integer,intent(out)::BoundaryList(2*(npx+npy)-4)
	real,intent(out)::BoundaryNodes(2*(npx+npy)-4,2)

	integer::InnerList(npInner),bottomList(npx),rightList(npy-1),upList(npx-1),leftList(npy-2)
	integer::i,n_inner,n_boundary,i_bottom,i_up,i_left,i_right
	real::dx,dy
    
    !initialization
	BoundaryList = 0
    n_inner = 0 ! inner domain nodes numbers
    n_boundary = 0 ! inner domain nodes numbers
	InnerList = 0
	bottomList = 0
	upList = 0
	leftList = 0
	rightList = 0
	i_bottom = 0
	i_up = 0
	i_left = 0
	i_right = 0
	dx = (x_end - x_start)/(npx-1)
	dy = (y_end - y_start)/(npy-1)
    
    write(*,*)"Boundary List generating..."    

    !*****************************************************************************************
    !***********  categorize the interior/boundary nodes
    !*****************************************************************************************
    do i=1,np 
        if ((xynodes(i,1) - x_start) > dx/10.  .and. (x_end - xynodes(i,1)) > dx/10. &
            .and. (xynodes(i,2) - y_start) > dy/10. .and. (y_end - xynodes(i,2)) > dy/10.) then
            n_inner = n_inner + 1       ! Inner domain nodes
            InnerList(n_inner) = i
        else
            n_boundary = n_boundary + 1 ! boundary nodes
            BoundaryList(n_boundary) = i
        end if        
    end do  
    
    if (n_inner .ne. npInner .or. n_boundary .ne. (np-npInner)) then
        write(*,*)"Error happens in whole Nodes Categorization. Exiting..."
        stop
    end if

    ! ****************************************************************************************
	!*********** This is to find the circular ordered EBC nodes list and its nodes coordinates
	! ****************************************************************************************
	do i=1,2*(npx+npy)-4
		if (xynodes(BoundaryList(i),2) - y_start < dy/10.) then  ! bottomList
			i_bottom = i_bottom + 1
			bottomList(i_bottom) = BoundaryList(i)
		
		else if(x_end - xynodes(BoundaryList(i),1) < dx/10.) then  ! rightList
			i_right = i_right + 1
			rightList(i_right) = BoundaryList(i)
			
		else if(y_end - xynodes(BoundaryList(i),2) < dy/10.) then  ! upList
			i_up = i_up + 1
			upList(i_up) = BoundaryList(i)
			
		else if(xynodes(BoundaryList(i),1) - x_start < dx/10.) then  ! leftList
			i_left = i_left + 1
			leftList(i_left) = BoundaryList(i)		
			
		else 
			write(*,*)"Error nodes on the boundary..."
			stop
			
		end if
		
	end do 
	
	if ((i_left .ne. (npy-2)) .or. (i_right .ne. (npy-1)) .or. (i_bottom .ne. npx) .or. (i_up .ne. (npx-1))) then
		write(*,*)"Error happens in boundary Nodes Categorization. Exiting..."
		stop
	end if
	
	BoundaryList(1:npx) =  bottomList
	BoundaryList(1+npx:npx+npy-1) =  rightList
	BoundaryList(npx+npy:2*npx+npy-2) = upList(npx-1:1:-1) ! up boundary: REVERSE!!!
	BoundaryList(2*npx+npy-1:2*(npx+npy)-4) = leftList(npy-2:1:-1) ! left boundary: REVERSE!!!
	
	BoundaryNodes = xynodes(BoundaryList(:),:)
	
    write(*,*)"Done!"
    
end subroutine

!********************************************************************************************************************

subroutine getOrderedNodeList(npx,npy,npEBC,EBCList,OrderedNodeList)
    implicit none
	! THIS IS FOR TRANSFORMATION METHOD USAGE
    ! NodeList = [ElseList,EBCNodeList]: from [1,2,...,np] ==> [...,EBCList]
	
	integer,intent(in)::npx,npy,npEBC
	integer,intent(in)::EBCList(npEBC)
	integer,intent(out)::OrderedNodeList(npx*npy)
	
	integer::NodeList(npx*npy),ElseList(npx*npy-npEBC),ElseListOrdered(npx*npy-npEBC)
	integer::EBCListOrdered(npEBC)
	integer::np,counter
	integer::i,j
	
    !!Initialization
	np = npx*npy
	OrderedNodeList = 0
	ElseList = 0
	EBCListOrdered = 0
	ElseListOrdered = 0
	counter = 0
	do i =1,np
		NodeList(i) = i
	end do
	
	! assign ElseList
	do i = 1, np
		do j = 1,npEBC 
			if (NodeList(i) == EBCList(j)) NodeList(i) = 0
		end do
	end do
	
	do i = 1,np
		if ( NodeList(i) /= 0) then
			counter = counter + 1
			ElseList(counter) = NodeList(i)
		end if
	end do
	
	! judge
	if (counter /= np-npEBC) stop "SOMETHING WRONG HAPPENS IN SUBROUTINE: getOrderedNodeList"
	
	call IncreaseOrderedList(npEBC,EBCList,EBCListOrdered)
	call IncreaseOrderedList(np-npEBC,ElseList,ElseListOrdered)
	
	OrderedNodeList(1:np-npEBC) = ElseListOrdered
	OrderedNodeList(1+np-npEBC:np) = EBCListOrdered

end subroutine

!********************************************************************************************************************

subroutine IncreaseOrderedList(n,List,ListOrdered)
! SORTING ALGORITHM: BUBBLE
	implicit none
	
	integer::n,List(n)
	integer::ListOrdered(n)
	
	integer::i,j
	
	ListOrdered = 0
	
	do i = 1,n
		do j = n,i+1,-1
			if (List(j-1) > List(j)) then ! swap
				List(j-1) = List(j-1) + List(j)
				List(j) = List(j-1) - List(j)
				List(j-1) = List(j-1) - List(j)
			end if
		end do
	end do
	
	ListOrdered = List
	
end subroutine

!********************************************************************************************************************



!********************************************************************************************************************