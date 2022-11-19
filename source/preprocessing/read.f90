subroutine read_control(x_start,x_end,y_start,y_end,npx,npy,nGx,nGy,nGBC,dilation,Emod,nu,Pforce,i_model,i_algorithm,beta)
	implicit none
	
	!###### problem domain ######		
	real::x_start, x_end
	real::y_start, y_end
	integer::npx,npy	 						! node numbers
	integer::nGx,nGy						 	! nGx*nGy Gauss points in a cell
	integer::nGBC								! n-points Gauss quadrature points in the Essential boundary	
	real::dilation,Emod,nu,beta,Pforce
	integer::i_model,i_algorithm
	
	integer::IERROR
	character(100)::ctemp
	
	open(222,file="control.dat",status="old")    

	do 
		read(222,'(A50)',iostat=IERROR)ctemp
		if(IERROR == -1) exit !EOF
		ctemp = trim(ctemp)
		
		if(ctemp == '*d_Xrange') then
			read(222,*) !comment
			read(222,*,iostat=IERROR)x_start,x_end
			if(IERROR /= 0) then
				stop "Error in reading d_Xrange. Exit..."
			end if
		end if
		
		if(ctemp == '*d_Yrange') then
			read(222,*) !comment
			read(222,*,iostat=IERROR)y_start,y_end
			if(IERROR /= 0) then
				stop "Error in reading d_Yrange. Exit..."
			end if
		end if
		
		if(ctemp == '*i_Nxynodes') then
			read(222,*) !comment
			read(222,*,iostat=IERROR)npx,npy
			if(IERROR /= 0) then
				stop "Error in reading i_Nxynodes. Exit..."
			end if
		end if
		
		if(ctemp == '*i_NxyGaussPts') then
			read(222,*) !comment
			read(222,*,iostat=IERROR)nGx,nGy
			if(IERROR /= 0) then
				stop "Error in reading i_NxyGaussPts. Exit..."
			end if
		end if
		
		if(ctemp == '*i_NGBC') then
			read(222,*) !comment
			read(222,*,iostat=IERROR)nGBC
			if(IERROR /= 0) then
				stop "Error in reading i_NGBC. Exit..."
			end if
		end if
		
		if(ctemp == '*d_Dilation') then
			read(222,*) !comment
			read(222,*,iostat=IERROR)dilation
			if(IERROR /= 0) then
				stop "Error in reading d_Dilation. Exit..."
			end if
		end if
		
		if(ctemp == '*d_Elastic') then
			read(222,*) !comment
			read(222,"(E8.2)",iostat=IERROR)Emod
			if(IERROR /= 0) then
				stop "Error in reading d_Elastic. Exit..."
			end if
		end if
		
		if(ctemp == '*d_Poisson') then
			read(222,*) !comment
			read(222,*,iostat=IERROR)nu
			if(IERROR /= 0) then
				stop "Error in reading d_Poisson. Exit..."
			end if
		end if
		
		if(ctemp == '*i_Model') then
			read(222,*) !comment
			read(222,*,iostat=IERROR)i_model
			if(IERROR /= 0) then
				stop "Error in reading i_Model. Exit..."
			end if
		end if
		
		if(ctemp == '*i_Algorithm') then
			read(222,*) !comment
			read(222,*,iostat=IERROR)i_algorithm
			if(IERROR /= 0) then
				stop "Error in reading i_algorithm. Exit..."
			end if
		end if
		
		if(ctemp == '*Pforce') then
			read(222,*) !comment
			read(222,"(E8.2)",iostat=IERROR)Pforce
			if(IERROR /= 0) then
				stop "Error in reading Pforce. Exit..."
			end if
		end if
		
		if(ctemp == '*d_Beta') then
			read(222,*) !comment
			read(222,*,iostat=IERROR)beta
			if(IERROR /= 0) then
				stop "Error in reading d_Beta. Exit..."
			end if
		end if
	
	end do

    close(222)	
	
end subroutine