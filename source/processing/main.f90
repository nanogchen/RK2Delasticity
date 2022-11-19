program RK2Delasticity_processing
	use constants
	implicit none
	
	real::xynodes(np,2),xyGaussPts(npGauss,2),xyGaussWts(npGauss)
	real::EBCnodes(npEBC,2),NBCnodes(npNBC,2)
	real::EBCGaussPts(npGaussEBC,2),NBCGaussPts(npGaussNBC,2)
	real::EBCGaussWts(npGaussEBC),NBCGaussWts(npGaussNBC)	
	integer::NeighborListInner(np,NumNeighborInner),NeighborListEBC(np,NumNeighborEBC),NeighborListNBC(np,NumNeighborNBC)
	
	real::mat_KBefo(2*np,2*np),Vec_fBefo(2*np),Vec_UBefo(2*np)	! common term
	real::mat_KAfter(2*np,2*np),Vec_fAfter(2*np),Vec_UAfter(2*np)
	real::InvMatrix_N(2*np,2*np)
	integer::OrderedNodeList(np)

	real::tipload(2) 
	
	integer::Trans_matRight(2*np,2*np)			! for transformation method
	real::Mat_A(2*npEBC,2*np),Vec_q(2*npEBC)	! Lagrangian multiplier
	real::solU(2*np),solLambda(2*npEBC)			! Lagrangian multiplier
	real::K_beta(2*np,2*np),f_beta(2*np) 		! for penalty method
	real::K_Nitsche(2*np,2*np),f_Nitsche(2*np) 	! for Nitsche's method
	
	real::uGaussApprox(2*npGauss),uGaussExact(2*npGauss)
	real::L2errorU,L2errorV,L2errorstress_xx
	real::strainGauss(3,npGauss), stressGauss_approx(3,npGauss), stressGauss_exact(3,npGauss)	
	
	real::time_begin,time_end
	
! !*******************************************************************************
! ! all the calculations, like FindKf and L2Norm, are done at Gauss Points
! !*******************************************************************************

	call CPU_TIME(time_begin)
	
	!Initialization
	xynodes = dble(0)
	xyGaussPts = dble(0)
	xyGaussWts = dble(0)
	EBCnodes = dble(0)
	EBCGaussPts = dble(0)
	EBCGaussWts = dble(0)
	NBCnodes = dble(0)
	NBCGaussPts = dble(0)
	NBCGaussWts = dble(0)
	OrderedNodeList = 0
	
	NeighborListInner = 0
	NeighborListEBC = 0
	NeighborListNBC = 0
	Trans_matRight = 0

	tipload = (/Pforce,0.0/)
	
	mat_KBefo = dble(0)
	Vec_fBefo = dble(0)
	mat_KAfter = dble(0)
	Vec_fAfter = dble(0)
	Vec_UBefo = dble(0)
	Vec_UAfter = dble(0)
	InvMatrix_N = dble(0)
	
	K_beta = dble(0)
	f_beta = dble(0)
	
	K_Nitsche = dble(0)
	f_Nitsche = dble(0)
	
	Mat_A = dble(0)
	Vec_q = dble(0)
	solU = dble(0)
	solLambda = dble(0)
	
	uGaussApprox = dble(0)
	uGaussExact = dble(0)	
	L2errorU = dble(0)
	L2errorV = dble(0)
	L2errorstress_xx = dble(0)

	strainGauss = dble(0)
	stressGauss_approx = dble(0)
	stressGauss_exact = dble(0)

!*************************************** calculating *****************************************
!*************************************** calculating *****************************************	
	call readInputData(xynodes,xyGaussPts,xyGaussWts,EBCnodes,NBCnodes,EBCGaussPts,NBCGaussPts,EBCGaussWts,NBCGaussWts,OrderedNodeList)
	call SearchNeighborList(xynodes,xyGaussPts,EBCGaussPts,NBCGaussPts,NeighborListInner,NeighborListEBC,NeighborListNBC)
	call FindKf(xynodes,xyGaussPts,xyGaussWts,NBCGaussPts,NBCGaussWts,NeighborListInner,NeighborListNBC,mat_KBefo,Vec_fBefo)
	!call AddTipLoad(np, tipload, Vec_fBefo) ! add the tip load at np-th node which is natural BC
	!call AddTipLoad(npx, tipload, Vec_fBefo) ! add the tip load at np-th node which is natural BC
	!call saveRealData2D(101, "mat_KBefoTrans.txt", 2*np, 2*np, mat_KBefo)

	if (i_algorithm == 1) then 		! transformation
		call transforMatrix(np,xynodes,InvMatrix_N,mat_KBefo)
		call CalExactDomainDisp(np,Vec_UBefo, xynodes) ! introduce the essential boundary condition
		call rearrangeMat(np,OrderedNodeList,mat_KBefo,Vec_fBefo,Vec_UBefo,mat_KAfter,Vec_fAfter,Vec_UAfter,Trans_matRight)
		call SolSysEq_transformation(np,npEBC,np-npEBC,Trans_matRight,mat_KAfter,Vec_fAfter,Vec_UAfter,solU)		
	
	else if	(i_algorithm == 2) then	 ! Lagrangian multiplier
		call FindAq(NeighborListEBC,xynodes,EBCnodes,EBCGaussPts,EBCGaussWts,Mat_A,Vec_q)
		call SolSysEq_LagMulti(mat_KBefo, Mat_A, Vec_fBefo, Vec_q, solU, solLambda)
		
	else if (i_algorithm == 3) then ! penalty method
		call FindKfBeta(xynodes,EBCGaussPts,EBCGaussWts,K_beta,f_beta)
		call SolSysEq_PenaltyMethod(np,beta,mat_KBefo,Vec_fBefo,K_beta,f_beta,solU)
		
	else if (i_algorithm == 4) then ! Nitsche's method
		call FindKfBeta(xynodes,EBCGaussPts,EBCGaussWts,K_beta,f_beta)
		call FindKfNitsche(xynodes,EBCGaussPts,EBCGaussWts,K_Nitsche,f_Nitsche)
		call SolSysEq_NitscheMethod(np,beta,mat_KBefo,Vec_fBefo,K_beta,f_beta,K_Nitsche,f_Nitsche,solU)

	else 
		stop "Wrong in i_algorithm!"
	end if

!************************************** strain at Gauss pts ***************************************************
!************************************** strain at Gauss pts ***************************************************
	call getStrain(np, npGauss, xynodes, xyGaussPts, solU, InvMatrix_N,strainGauss)
	call saveRealData1D(121, "strain_xx_a.txt", npGauss, strainGauss(1,:))

!************************************** stress at Gauss pts ***************************************************
!************************************** stress at Gauss pts ***************************************************
	call getStress(npGauss, strainGauss,stressGauss_approx)
	call saveRealData1D(122, "stress_xx_a.txt", npGauss, stressGauss_approx(1,:))

	call CalExactDomainStress(npGauss, xyGaussPts, stressGauss_exact)
	call saveRealData1D(123, "stress_xx_e.txt", npGauss, stressGauss_exact(1,:))	

	call L2Norm(L2errorstress_xx,npGauss,xyGaussWts,stressGauss_exact(1,:),stressGauss_approx(1,:))
	call saveRealData1D(211, "L2stress_xx.txt", 1, L2errorstress_xx)

	
!************************************** L2 norm ***************************************************
!************************************** L2 norm ***************************************************
	call CalApproxValue_transformation(np,xynodes,InvMatrix_N,solU,npGauss,xyGaussPts,uGaussApprox)
	call CalExactDomainDisp(npGauss,uGaussExact, xyGaussPts)
	call L2Norm(L2errorU,npGauss,xyGaussWts,uGaussExact(1:2*npGauss:2),uGaussApprox(1:2*npGauss:2))
	call L2Norm(L2errorV,npGauss,xyGaussWts,uGaussExact(2:2*npGauss:2),uGaussApprox(2:2*npGauss:2))
	
!*************************************** for post-processing *****************************************
!*************************************** for post-processing *****************************************	
	call constantFileWrite(x_start,x_end,y_start,y_end,npx,npy,nGx,nGy,nGBC,npEBC,npNBC,npGaussEBC,npGaussNBC, &
							dilation,Emod,nu,i_algorithm,beta)
	
!*************************************** saving the data *****************************************
!*************************************** saving the data *****************************************
	call saveRealData2D(11, "mat_KBefo.txt", 2*np, 2*np, mat_KBefo)
	call saveRealData1D(12, "Vec_fBefo.txt", 2*np, Vec_fBefo)
	call saveRealData1D(13, "Vec_uBefo(exact).txt", 2*np, Vec_uBefo)
	call saveRealData1D(14,"solU.txt",2*np,solU)
	
	call vtkWriteDisp(np,xynodes,solU) ! Nodal coordinate and solution
	
	if (i_algorithm == 1) then	! transformation
		call saveRealData1D(141, "uGaussApprox.txt", 2*npGauss, uGaussApprox)
		call saveRealData1D(142, "uGaussExact.txt", 2*npGauss, uGaussExact)
		call saveRealData1D(143, "L2errorU.txt", 1, L2errorU)
		call saveRealData1D(144, "L2errorV.txt", 1, L2errorV)
	
	else if (i_algorithm == 2) then	! Lagrangian multiplier
		call saveRealData2D(15, "mat_A.txt", 2*npEBC,2*np, mat_A)
		call saveRealData1D(16, "Vec_q.txt", 2*npEBC, Vec_q)	
		call saveRealData1D(17,"solLambda.txt",2*npEBC,solLambda)
	
	else if (i_algorithm == 3) then	! penalty method
		call saveRealData2D(18, "K_beta.txt", 2*np,2*np, K_beta)
		call saveRealData1D(19, "f_beta.txt", 2*np, f_beta)
	
	else if (i_algorithm == 4) then	! Nitsche's method
		call saveRealData2D(20, "K_beta.txt", 2*np,2*np, K_beta)
		call saveRealData1D(21, "f_beta.txt", 2*np, f_beta)
		call saveRealData2D(22, "K_Nitsche.txt", 2*np,2*np, K_Nitsche)
		call saveRealData1D(23, "f_Nitsche.txt", 2*np, f_Nitsche)
	
	else
		stop "Wrong in i_algorithm!"
	end if	
	
	
	call CPU_TIME(time_end)
	write(*,*) "Computation time is:", time_end-time_begin
	
end program
