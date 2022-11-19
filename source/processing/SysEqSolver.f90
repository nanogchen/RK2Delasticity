!*******************************************************************************

subroutine SolSysEq_transformation(np,npEBC,npElse,Trans_matRight,matAfter,FAfter,Vec_UAfter,solU) 
! this is to solve the system equation, contains the procedure to make block matrix and solve it
	implicit none
	
	integer,intent(in)::np,npEBC,npElse,Trans_matRight(2*np,2*np)
	real,intent(in)::matAfter(2*np,2*np),FAfter(2*np),Vec_UAfter(2*np)
	real,intent(inout)::solU(2*np) ! this is the solution vector of nodal value: solU[1,np]
	
	real::Kaa(2*npElse,2*npElse),Kab(2*npElse,2*npEBC),Kba(2*npEBC,2*npElse),Kbb(2*npEBC,2*npEBC)
	real::InvKaa(2*npElse,2*npElse),f_temp(2*npElse)
	real::fa(2*npElse),ua(2*npElse),ub(2*npEBC)
	
	!initialization
	Kaa = dble(0)
	Kab = dble(0)
	Kba = dble(0)
	Kbb = dble(0)
	fa = dble(0) 
	! fb = dble(0) ! of no use
	f_temp = dble(0)
	InvKaa = dble(0)
	
	write(*,*)"Solving the ultimate system equation by transformation method..."
	
	! [Kaa,Kab;Kba,Kbb]*[Ua;Ub]=[fa;fb]
	
	fa = FAfter(1:2*npElse)
	ua = Vec_UAfter(1:2*npElse)
	ub = Vec_UAfter(2*npElse+1:2*np)
	
	Kaa = matAfter(1:2*npElse,1:2*npElse)
	Kab = matAfter(1:2*npElse,2*npElse+1:2*np)
	Kba = matAfter(2*npElse+1:2*np, 1:2*npElse)
	Kbb = matAfter(2*npElse+1:2*np, 2*npElse+1:2*np)
	
	call inverse(Kaa,2*npElse,InvKaa)
	f_temp = fa - matmul(Kab,ub)
	ua = matmul(InvKaa,f_temp)
	
	solU(1:2*npElse) = ua
	solU(2*npElse+1:2*np) = ub
	
	solU = matmul(Trans_matRight,solU)
	
	write(*,*)"Done!"
	
end subroutine

!*******************************************************************************

subroutine SolSysEq_NitscheMethod(np,beta,mat_k,vec_f,K_beta,f_beta,K_Nitsche,f_Nitsche,solU)
	!	This is to assemble and solve the system equation obtained by Nitsche method.
	implicit none
	
	integer,intent(in)::np
	real(kind=16),intent(in)::beta
	real,intent(in)::mat_k(2*np,2*np),vec_f(2*np)   ! mat_k,vec_f are calculated by FindKf
	real,intent(in)::K_beta(2*np,2*np),f_beta(2*np) ! K_beta,f_beta are calculated by FindKfBeta
	real,intent(in)::K_Nitsche(2*np,2*np),f_Nitsche(2*np) ! K_Nitsche,f_Nitsche are calculated by FindKfNitsche
	real,intent(out)::solU(2*np)
	
	real::KK(2*np,2*np),bb(2*np) ! temp variables
	real::invKK(2*np,2*np)
	
	write(*,*)"Solving the ultimate system equation by Nitsche's method..."
	
	!Initialization
	solU = dble(0)
	KK = dble(0)
	bb = dble(0)
	invKK = dble(0)
	
	! assemble the system Eq.: KK*U=bb
	KK = mat_k + beta*K_beta - K_Nitsche - transpose(K_Nitsche)
	bb = vec_f + beta*f_beta - f_Nitsche
	
	! find the solution
	call inverse(KK,2*np,invKK)
	solU = matmul(invKK,bb)
	
	write(*,*)"Done!"
	
end subroutine

!*******************************************************************************

subroutine SolSysEq_PenaltyMethod(np,beta,mat_k,vec_f,K_beta,f_beta,solU)
!	This is to assemble and solve the system equation obtained by penalty method.
	implicit none
	
	integer,intent(in)::np
	real(kind=16),intent(in)::beta
	real,intent(in)::mat_k(2*np,2*np),vec_f(2*np),K_beta(2*np,2*np),f_beta(2*np) ! mat_k,vec_f are calculated by FindKf
	real,intent(out)::solU(2*np)
	
	real::KK(2*np,2*np),bb(2*np) ! temp variables
	real::invKK(2*np,2*np)
	
	write(*,*)"Solving the ultimate system equation by Penalty method..."
	
	!Initialization
	solU = dble(0)
	KK = dble(0)
	bb = dble(0)
	invKK = dble(0)
	
	! assemble the system Eq.: KK*U=bb
	KK = mat_k + beta*K_beta
	bb = vec_f + beta*f_beta
	
	! find the solution
	call inverse(KK,2*np,invKK)
	solU = matmul(invKK,bb)
	
	write(*,*)"Done!"
	
end subroutine

!*******************************************************************************

subroutine SolSysEq_LagMulti(Mat_K, Mat_A, Vec_f, Vec_q, solU, solLambda)
!	This is to assemble and solve the system equation obtained by Lagrangian Multiplier.
	use constants,only:np, npEBC
	implicit none
	
	real,intent(in)::Mat_K(2*np,2*np), Vec_f(2*np), Mat_A(2*npEBC,2*np),Vec_q(2*npEBC) ! Mat_K,Vec_f are calculated by FindKf
	real,intent(out)::solU(2*np), solLambda(2*npEBC)
	
	real::tempK(2*(np+npEBC),2*(np+npEBC)),tempX(2*(np+npEBC)),tempB(2*(np+npEBC)) ! assemble the KX=B sys Eq.
	real::invK(2*(np+npEBC),2*(np+npEBC))
	real::Trans_Mat_A(2*np,2*npEBC)
	
	write(*,*)"Solving the ultimate system equation by Lagrangian Multiplier..."
	
	!Initialization
	solU = dble(0)
	solLambda = dble(0)
	tempK = dble(0)
	invK = dble(0)
	tempX = dble(0)
	tempB = dble(0)
	Trans_Mat_A = transpose(Mat_A)
	
	! assemble sys Eq.
	tempK(1:2*np,1:2*np) = Mat_K
	tempK(1:2*np,2*np+1:2*(np+npEBC)) = Trans_Mat_A
	tempK(2*np+1:2*(np+npEBC),1:2*np) = Mat_A
	
	tempB(1:2*np) = Vec_f
	tempB(2*np+1:2*(np+npEBC)) = Vec_q
	
	! solving Eq.
	call inverse(tempK,2*(np+npEBC),invK)
	tempX = matmul(invK,tempB)
	
	solU = tempX(1:2*np)
	solLambda = tempX(2*np+1:2*(np+npEBC))	
	
	write(*,*)"Done!"
	
end subroutine

!*******************************************************************************