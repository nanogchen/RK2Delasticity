subroutine SplitDisp(np,solU,dispU,dispV)
	implicit none
	
	integer,intent(in)::np
	real,intent(in)::solU(2*np)
	real,intent(out)::dispU(np),dispV(np)
	
	integer::i
	
!	Initialization
	dispU = dble(0)
	dispV = dble(0)
	
	do i=1,np
		dispU(i) = solU(2*(i-1)+1)
		dispV(i) = solU(2*(i-1)+2)
	end do 
	
end subroutine