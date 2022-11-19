!*******************************************************************************

subroutine FindKf(xynodes,xyGaussPts,xyGaussWts,NBCGaussPts,NBCGaussWts,NeighborListInner,NeighborListNBC,mat_K,Vec_f)
!	this is to find the direct matrix K and f for the system equation Ku=f
	use constants,only:np,npGauss,npNBC,npGaussNBC,NumNeighborInner,NumNeighborNBC,i_model,HookMat_strain, HookMat_stress
	implicit none
	
	real,intent(in)::xynodes(np,2),xyGaussPts(npGauss,2),xyGaussWts(npGauss)
	real,intent(in)::NBCGaussPts(npGaussNBC,2),NBCGaussWts(npGaussNBC)
	integer,intent(in)::NeighborListInner(np,NumNeighborInner),NeighborListNBC(np,NumNeighborNBC)
	real,intent(out)::mat_K(2*np,2*np),Vec_f(2*np)
	
	real::HookMat(3,3)
	real::Ni, Ndxyi(2), Nj, Ndxyj(2), Ui, force(2), Nf, Ndxyf(2),Kij(2,2)
	real::Bi(3,2),Bj(3,2),PhiI(2,2),PhiJ(2,2)
	real::N_mat(2,2)
	integer::i,j,k,counter	
	
	!initialization
	Bi = dble(0)
	Bj = dble(0)
	PhiI = dble(0)
	PhiJ = dble(0)
	mat_K = dble(0)
	Kij = dble(0)
	Vec_f = dble(0)
	Ni = dble(0)
	Ndxyi = dble(0)
	Nj = dble(0)
	Ndxyj = dble(0)
	force = dble(0)
	Ui = dble(0)  ! for nodal value of essential nodes
	Nf = dble(0)
	Ndxyf = dble(0)
	N_mat = dble(0)
	HookMat = dble(0)
	
	write(*,*)"Finding K, f for the system equation..."
	
	if (i_model == 1) then
		HookMat = HookMat_stress
	else if (i_model == 2) then
		HookMat = HookMat_strain
	else 
		stop "Something wrong happen in choosing the plane stress/strain model."
	end if
	
	do i=1,np
		do j=1,np  
			do counter=1, NumNeighborInner
				if (NeighborListInner(i,counter) .ne. 0) then
					k = NeighborListInner(i,counter)
					call CalN(Ni, Ndxyi, xynodes, xynodes(i,1), xynodes(i,2), xyGaussPts(k,1), xyGaussPts(k,2))
					call CalN(Nj, Ndxyj, xynodes, xynodes(j,1), xynodes(j,2), xyGaussPts(k,1), xyGaussPts(k,2))
					
					Bi = reshape((/ Ndxyi(1),0.0, Ndxyi(2),0.0,Ndxyi(2),Ndxyi(1)/),(/3,2/))
					Bj = reshape((/ Ndxyj(1),0.0, Ndxyj(2),0.0,Ndxyj(2),Ndxyj(1)/),(/3,2/))
					
					Kij = matmul(matmul(transpose(Bi),HookMat),Bj)*xyGaussWts(k) ! mat_K(i,j) is 2*2 matrix 
					
					mat_K((i-1)*2+1,(j-1)*2+1) = mat_K((i-1)*2+1,(j-1)*2+1) + Kij(1,1) ! mat_K(i,j) is 2*2 matrix
					mat_K((i-1)*2+1,(j-1)*2+2) = mat_K((i-1)*2+1,(j-1)*2+2) + Kij(1,2)
					mat_K(i*2,(j-1)*2+1) = mat_K(i*2,(j-1)*2+1) + Kij(2,1)
					mat_K(i*2,j*2) = mat_K(i*2,j*2) + Kij(2,2)
					
				end if
			end do
		end do 
	end do
	
! ! on the domain and the natural BC
	! ! on the domain
	do i=1,np
		do counter=1,NumNeighborInner
			k = NeighborListInner(i,counter) 
			if (k .ne. 0) then
				call CalN(Nf, Ndxyf, xynodes, xynodes(i,1), xynodes(i,2), xyGaussPts(k,1), xyGaussPts(k,2))
				call CalBodyForce(force,xyGaussPts(k,:))
				
				N_mat = reshape((/ Nf,0.0,0.0,Nf /),(/2,2/))
				Vec_f(2*(i-1)+1:2*(i-1)+2) = Vec_f(2*(i-1)+1:2*(i-1)+2) + matmul(N_mat,force)*xyGaussWts(k)
			end if
		end do
	end do
	
	! on the natural BC
	if (npNBC > 0) then
		do i=1,np
			do counter=1,NumNeighborNBC
				k = NeighborListNBC(i,counter)
				if (k .ne. 0) then 
					call CalN(Nf, Ndxyf, xynodes, xynodes(i,1), xynodes(i,2), NBCGaussPts(k,1), NBCGaussPts(k,2))
					
					PhiI = reshape((/ Nf, 0., 0., Nf/),(/2,2/))					
					call CalTraction(force,NBCGaussPts(k,:))
					
					Vec_f(2*(i-1)+1:2*(i-1)+2) = Vec_f(2*(i-1)+1:2*(i-1)+2) + matmul(PhiI,force)*NBCGaussWts(k)
				end if
			end do
		end do
	end if
	
	write(*,*)"Done!"
	
end subroutine

!*******************************************************************************

subroutine FindAq(NeighborListEBC,xynodes,EBCnodes,EBCGaussPts,EBCGaussWts,Mat_A,Vec_q)
!	This is to find other matrics and vectors in system equation by Lagrangian Multiplier.
	use constants,only:npEBC,npGaussEBC,np,NumNeighborEBC
	implicit none
	
	integer,intent(in)::NeighborListEBC(np,NumNeighborEBC)
	real,intent(in)::xynodes(np,2),EBCnodes(npEBC,2),EBCGaussPts(npGaussEBC,2),EBCGaussWts(npGaussEBC)
	real,intent(out)::Mat_A(2*npEBC,2*np),Vec_q(2*npEBC)
	
	real::Nfe,N,Ndxy(2),disp_u(2)
	integer::i,j,k,counter
	
	! matrix form
	real::Nfe_mat(2,2),N_mat(2,2),Mat_Aij(2,2),Vec_qi(2)
	
	write(*,*)"Finding A, q for the system equation..."
	
	!Initialization
	Mat_A = dble(0) 
	Vec_q = dble(0) 
	
	Nfe = dble(0)
	N = dble(0)
	Ndxy = dble(0)
	disp_u = dble(0)
	
	Nfe_mat = dble(0)
	N_mat = dble(0)
	Mat_Aij = dble(0)
	Vec_qi = dble(0)
	
	if (npEBC > 0) then
		do i=1,npEBC
			do j=1,np
				do counter=1,NumNeighborEBC
					k = NeighborListEBC(i,counter)
					if (k .ne. 0) then
						call CalN_FE(Nfe,EBCnodes(i,:),EBCGaussPts(k,:))
						call CalN(N, Ndxy, xynodes, xynodes(j,1), xynodes(j,2), EBCGaussPts(k,1), EBCGaussPts(k,2))
						
						Nfe_mat = reshape((/ Nfe,0.0,0.0,Nfe /),(/2,2/))
						N_mat = reshape((/ N,0.0,0.0,N /),(/2,2/))
						
						Mat_Aij = matmul(Nfe_mat,transpose(N_mat))*EBCGaussWts(k)
						
						Mat_A((i-1)*2+1,(j-1)*2+1) = Mat_A((i-1)*2+1,(j-1)*2+1) + Mat_Aij(1,1)
						Mat_A((i-1)*2+1,(j-1)*2+2) = Mat_A((i-1)*2+1,(j-1)*2+2) + Mat_Aij(1,2)
						Mat_A(i*2,(j-1)*2+1) = Mat_A(i*2,(j-1)*2+1) + Mat_Aij(2,1)
						Mat_A(i*2,j*2) = Mat_A(i*2,j*2) + Mat_Aij(2,2)
						
						! Mat_A((i-1)*2+1:(i-1)*2+2,(j-1)*2+1:(j-1)*2+2) = Mat_A((i-1)*2+1:(i-1)*2+2,(j-1)*2+1:(j-1)*2+2) + Mat_Aij
					endif
				end do
			end do 
		end do
		
		do i = 1,npEBC
			do counter=1,NumNeighborEBC
				k = NeighborListEBC(i,counter)
				if (k .ne. 0) then			
					call CalN_FE(Nfe,EBCnodes(i,:),EBCGaussPts(k,:))			
					Nfe_mat = reshape((/ Nfe,0.0,0.0,Nfe /),(/2,2/))
					
					call CalExactDisp(disp_u,EBCGaussPts(k,:))
					
					Vec_qi = matmul(Nfe_mat,disp_u)*EBCGaussWts(k)			
					Vec_q(2*(i-1)+1:2*(i-1)+2) = Vec_q(2*(i-1)+1:2*(i-1)+2) + Vec_qi
				endif
			end do
		end do	
		
	end if	
	write(*,*)"Done!"	
	
end subroutine

!*******************************************************************************

subroutine FindKfBeta(xynodes,EBCGaussPts,EBCGaussWts,K_beta,f_beta)
	! This is to find the sys matrix and vector in penalty method
	use constants,only:np,npGauss,npGaussEBC
	implicit none
	
	real,intent(in)::xynodes(np,2),EBCGaussPts(npGaussEBC,2),EBCGaussWts(npGaussEBC)
	real,intent(out)::K_beta(2*np,2*np),f_beta(2*np)
	
	integer::i,j,k
	real::N,Ndxy(2),Ni,Ndxyi(2),Nj,Ndxyj(2),disp_u(2)
	real::mat_Ni(2,2),mat_Nj(2,2),mat_N(2,2)
	
	write(*,*)"Finding K_beta, f_beta for the system equation..."
	
	!Initialization
	K_beta = dble(0)
	f_beta = dble(0)
	N = dble(0)
	Ndxy = dble(0)
	Ni = dble(0)
	Nj = dble(0)
	Ndxyi = dble(0)
	Ndxyj = dble(0)
	disp_u = dble(0)
	mat_Ni = dble(0)
	mat_Nj = dble(0)
	mat_N = dble(0)
	
	do i=1,np 
		do j=1,np
			do k = 1,npGaussEBC
				call CalN(Ni, Ndxyi, xynodes, xynodes(i,1), xynodes(i,2), EBCGaussPts(k,1), EBCGaussPts(k,2))
				call CalN(Nj, Ndxyj, xynodes, xynodes(j,1), xynodes(j,2), EBCGaussPts(k,1), EBCGaussPts(k,2))
				
				mat_Ni = reshape((/Ni, 0.0, 0.0, Ni/),(/2,2/))
				mat_Nj = reshape((/Nj, 0.0, 0.0, Nj/),(/2,2/))
				
				K_beta(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2) = K_beta(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2) &
																	+ matmul(mat_Ni,mat_Nj)*EBCGaussWts(k)
			end do 
		end do 
		
		do k = 1, npGaussEBC
			call CalN(N, Ndxy, xynodes, xynodes(i,1), xynodes(i,2), EBCGaussPts(k,1), EBCGaussPts(k,2))
			call CalExactDisp(disp_u, EBCGaussPts(k,:))
			
			mat_N = reshape((/N, 0.0, 0.0, N/),(/2,2/))
			
			f_beta(2*(i-1)+1:2*(i-1)+2) = f_beta(2*(i-1)+1:2*(i-1)+2) + matmul(mat_N,disp_u)*EBCGaussWts(k)
		end do 
	
	end do
	
	write(*,*)"Done!"
	
end subroutine

!*******************************************************************************

subroutine FindKfNitsche(xynodes,EBCGaussPts,EBCGaussWts,K_Nitsche,f_Nitsche)
	use constants,only:np,npGaussEBC,x_start,x_end,y_start,y_end,dx,dy
	! This is to find the sys matrix and vector in Nitsche's method
	implicit none
	
	real,intent(in)::xynodes(np,2),EBCGaussPts(npGaussEBC,2),EBCGaussWts(npGaussEBC)
	real,intent(out)::K_Nitsche(2*np,2*np),f_Nitsche(2*np)
	
	integer::i,j,k
	real::N,Ndxy(2),Ni,Ndxyi(2),Nj,Ndxyj(2),disp_u(2)
	real::mat_Ni(2,2),mat_Nj(2,2),mat_N(2,2)
	real::norm(3,2),Bi(3,2),Bj(3,2),B(3,2)
	integer::unitnorm(2)
	
	write(*,*)"Finding K_Nitsche,f_Nitsche for the system equation..."
	
	!Initialization
	K_Nitsche = dble(0)
	f_Nitsche = dble(0)
	N = dble(0)
	Ndxy = dble(0)
	Ni = dble(0)
	Nj = dble(0)
	Ndxyi = dble(0)
	Ndxyj = dble(0)
	disp_u = dble(0)
	mat_Ni = dble(0)
	mat_Nj = dble(0)
	mat_N = dble(0)
	norm = dble(0)
	Bi = dble(0)
	Bj = dble(0)
	B = dble(0)
	
	do i=1,np 
		do j=1,np
			do k = 1,npGaussEBC
				call CalN(Ni, Ndxyi, xynodes, xynodes(i,1), xynodes(i,2), EBCGaussPts(k,1), EBCGaussPts(k,2))
				call CalN(Nj, Ndxyj, xynodes, xynodes(j,1), xynodes(j,2), EBCGaussPts(k,1), EBCGaussPts(k,2))
				
				mat_Ni = reshape((/Ni, 0.0, 0.0, Ni/),(/2,2/))
				Bj = reshape((/ Ndxyj(1),0.0, Ndxyj(2),0.0,Ndxyj(2),Ndxyj(1)/),(/3,2/))								
				
				
				if (abs(EBCGaussPts(k,2) - y_start) < dy/100.) then ! n = (0,-1)
					unitnorm = (/0,-1/) !  n = (nx,ny)
					norm = reshape((/unitnorm(1),0,unitnorm(2),0,unitnorm(2),unitnorm(1)/),(/3,2/))
					K_Nitsche(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2) = K_Nitsche(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2) &
												+ matmul(mat_Ni,matmul(transpose(Bj),norm))*EBCGaussWts(k)
				end if 
				
				if (abs(EBCGaussPts(k,1) - x_end) < dx/100.) then ! n = (1,0)
					unitnorm = (/1,0/) !  n = (nx,ny)
					norm = reshape((/unitnorm(1),0,unitnorm(2),0,unitnorm(2),unitnorm(1)/),(/3,2/))
					K_Nitsche(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2) = K_Nitsche(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2)&
												+ matmul(mat_Ni,matmul(transpose(Bj),norm))*EBCGaussWts(k)
				end if
				
				if (abs(EBCGaussPts(k,2) - y_end) < dy/100.) then ! n = (0,1)
					unitnorm = (/0,1/) !  n = (nx,ny)
					norm = reshape((/unitnorm(1),0,unitnorm(2),0,unitnorm(2),unitnorm(1)/),(/3,2/))
					K_Nitsche(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2) = K_Nitsche(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2)&
												+ matmul(mat_Ni,matmul(transpose(Bj),norm))*EBCGaussWts(k)
				end if 
				
				if (abs(EBCGaussPts(k,1) - x_start) < dx/100.) then ! n = (-1,0)
					unitnorm = (/-1,0/) !  n = (nx,ny)
					norm = reshape((/unitnorm(1),0,unitnorm(2),0,unitnorm(2),unitnorm(1)/),(/3,2/))
					K_Nitsche(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2) = K_Nitsche(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2)&
												+ matmul(mat_Ni,matmul(transpose(Bj),norm))*EBCGaussWts(k)
				end if

			end do 
		end do 
		
		do k = 1, npGaussEBC
			call CalN(N, Ndxy, xynodes, xynodes(i,1), xynodes(i,2), EBCGaussPts(k,1), EBCGaussPts(k,2))
			call CalExactDisp(disp_u, EBCGaussPts(k,:))

			B = reshape((/ Ndxy(1),0.0, Ndxy(2),0.0,Ndxy(2),Ndxy(1)/),(/3,2/))					
			
			if (abs(EBCGaussPts(k,2) - y_start) < dy/100.) then ! n = (0,-1)
				unitnorm = (/0,-1/) !  n = (nx,ny)
				norm = reshape((/unitnorm(1),0,unitnorm(2),0,unitnorm(2),unitnorm(1)/),(/3,2/))
				
				f_Nitsche(2*(i-1)+1:2*(i-1)+2) = f_Nitsche(2*(i-1)+1:2*(i-1)+2) + &
											matmul(matmul(transpose(B),norm),disp_u)*EBCGaussWts(k)
			end if 
			
			if (abs(EBCGaussPts(k,1) - x_end) < dx/100.) then ! n = (1,0)
				unitnorm = (/1,0/) !  n = (nx,ny)
				norm = reshape((/unitnorm(1),0,unitnorm(2),0,unitnorm(2),unitnorm(1)/),(/3,2/))
				
				f_Nitsche(2*(i-1)+1:2*(i-1)+2) = f_Nitsche(2*(i-1)+1:2*(i-1)+2) + &
											matmul(matmul(transpose(B),norm),disp_u)*EBCGaussWts(k)
			end if   						
			
			if (abs(EBCGaussPts(k,2) - y_end) < dy/100.) then ! n = (0,1)
				unitnorm = (/0,1/) !  n = (nx,ny)
				norm = reshape((/unitnorm(1),0,unitnorm(2),0,unitnorm(2),unitnorm(1)/),(/3,2/))
				
				f_Nitsche(2*(i-1)+1:2*(i-1)+2) = f_Nitsche(2*(i-1)+1:2*(i-1)+2) + &
											matmul(matmul(transpose(B),norm),disp_u)*EBCGaussWts(k)
			end if 
			
			if (abs(EBCGaussPts(k,1) - x_start) < dx/100.) then ! n = (-1,0)
				unitnorm = (/-1,0/) !  n = (nx,ny)
				norm = reshape((/unitnorm(1),0,unitnorm(2),0,unitnorm(2),unitnorm(1)/),(/3,2/))
				
				f_Nitsche(2*(i-1)+1:2*(i-1)+2) = f_Nitsche(2*(i-1)+1:2*(i-1)+2) + &
											matmul(matmul(transpose(B),norm),disp_u)*EBCGaussWts(k)
											
			end if
				
		end do 
	
	end do
	
end subroutine

!*******************************************************************************
