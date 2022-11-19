program RK2Delasticity_postprocessing
	use constants
	implicit none
	
	
	real::xynodes(np,2),xyGaussPts(npGauss,2),xyGaussWts(npGauss),solU(2*np),ExactsolU(2*np)
	real::dispU(np),dispV(np)
	real::exact_dispU(np),exact_dispV(np)
	real::xpts(npx),ypts(npy)
	integer::NodeNumbers(3)
	
	NodeNumbers = (/np,npx,npy/)
	
!*************************************** post-processing *****************************************
!*************************************** post-processing *****************************************		
	call readSolData(xynodes,xyGaussPts,xyGaussWts,solU,ExactsolU)
	call SplitDisp(np,solU,dispU,dispV)
	call SplitDisp(np,ExactsolU,exact_dispU,exact_dispV)
	
!*****************generate the coordinates for matlab ******************************
!*****************generate the coordinates for matlab ******************************
	call generate1DPts(x_start,x_end,npx,xpts)
	call generate1DPts(y_start,y_end,npy,ypts)
	
!*************************************** saving data *****************************************
!*************************************** saving data *****************************************
	call saveIntData1D(1, "NodeNumbers.txt", 3, NodeNumbers)
	call saveRealData1D(2, "xpts.txt", npx, xpts)
	call saveRealData1D(3, "ypts.txt", npy, ypts)
	call saveRealData1D(4, "dispU.txt", np, dispU)
	call saveRealData1D(5, "dispV.txt", np, dispV)
	call saveRealData1D(66, "exact_dispU.txt", np, exact_dispU)
	call saveRealData1D(7, "exact_dispV.txt", np, exact_dispV)	
	
	
	
end program