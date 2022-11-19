preprocessing: it is a general geometry generator, no need to change.

processing: only for elasticity, and the boundary condition should be altered adjusted to specific problem.
	e.g.: displacement at essential boundary, and traction force at natural boundary.
	it can be solved as input (in the form of array)
	
	read(1111,file="./data_pre/infile.txt")   ! xynodes, xyGaussPts, xyGaussWts, EBCnodes, NBCnodes ...
	write(2222,file="./data_pro/outfile.txt") ! systemMatrix, solU
	
  - processing_elastic: solver for 2D elastic problem
  - processing_pde: solver for 2D pde
  
postprocessing:
	there is no need to split the displacement vector, just in 2*n form. The only reason to do it is to 
	show the visualization in MATLAB.