subroutine constantFileWrite(x_start,x_end,y_start,y_end,npx,npy,nGx,nGy,nGBC,npEBC,npNBC,npGaussEBC,npGaussNBC,&
								dilation,Emod,nu,i_algorithm,beta)
	implicit none
	
	real,intent(in)::x_start,x_end,y_start,y_end,dilation,Emod,nu,beta
	integer,intent(in)::npx,npy,nGx,nGy,nGBC,npEBC,npNBC,npGaussEBC,npGaussNBC,i_algorithm
	
	open(10,file="../../source/postprocessing/constants.f95")
	write(10,*)"module constants"
	write(10,*)"	implicit none"	
	write(10,*)
	write(10,*)"	!domain geometry"
	write(10,"(A30,F4.1)")"	real,parameter::x_start = ",x_start
	write(10,"(A28,F4.1)")"	real,parameter::x_end = ",x_end
	write(10,"(A30,F4.1)")"	real,parameter::y_start = ",y_start
	write(10,"(A28,F4.1)")"	real,parameter::y_end = ",y_end
	write(10,*)"	real,parameter::Length = x_end - x_start"
	write(10,*)"	real,parameter::Height = y_end - y_start"	
	write(10,*)
	write(10,*)"	!domain node/Guass point generating parameter"
	write(10,"(A29,I3)")"integer,parameter::npx = ",npx
	write(10,"(A29,I3)")"integer,parameter::npy = ",npy
	write(10,*)"	integer,parameter::np = npx*npy"
	write(10,"(A29,I3)")"integer,parameter::nGx = ",nGx
	write(10,"(A29,I3)")"integer,parameter::nGy = ",nGy
	write(10,*)"	integer,parameter::npxGauss = nGx*(npx-1)"
	write(10,*)"	integer,parameter::npyGauss = nGy*(npy-1)"
	write(10,*)"	integer(kind=4),parameter::npGauss = npxGauss*npyGauss"	
	write(10,*)
	write(10,*)"	!boundary Guass point generating parameter"
	write(10,"(A30,I3)")"integer,parameter::nGBC = ", nGBC
	write(10,*)"	integer,parameter::npBC = 2*(npx+npy)-4"
	write(10,*)"	integer,parameter::npInner = np - npBC"
	write(10,"(A31,I3)")"integer,parameter::npEBC = ", npEBC
	write(10,"(A31,I3)")"integer,parameter::npNBC = ", npNBC
	write(10,"(A36,I3)")"integer,parameter::npGaussEBC = ", npGaussEBC
	write(10,"(A36,I3)")"integer,parameter::npGaussNBC = ", npGaussNBC
	write(10,*)
	write(10,*)"	!parameter for penalty method"
	write(10,"(A27,E8.2)")"real,parameter::beta = ", beta
	write(10,"(A37,I3)")"integer,parameter::i_algorithm = ", i_algorithm
	write(10,*)
	write(10,*)"	!Nodal configuration"
	write(10,"(A31,F4.2)")"real,parameter::dilation = ", dilation
	write(10,*)"	real,parameter::dx = (x_end - x_start)/(npx-1)"
	write(10,*)"	real,parameter::dy = (y_end - y_start)/(npy-1)"
	write(10,*)"	real,parameter::hx = dilation*dx ! for rectangular phi"
	write(10,*)"	real,parameter::hy = dilation*dy ! for rectangular phi"
	write(10,*)"	real,parameter::hr = dilation*sqrt(dy**2+dx**2) ! for circular phi"
	write(10,*)
	write(10,*)"	!Neighbor list to search neighbor Gauss points for each nodes"
	write(10,*)"	integer,parameter::NumNeighborInner = (( int(dilation) + 1 )*2)**2*nGx*nGy"
	write(10,*)"	integer,parameter::NumNeighborEBC = 8*nGBC ! the maximum num of neighbors for an essential boundary node"
	write(10,*)"	integer,parameter::NumNeighborNBC = 8*nGBC ! the maximum num of neighbors for a natural boundary node"
	write(10,*)
	write(10,*)"	!material properties"
	write(10,*)"	real,parameter::Imoment = Height**3/12.0"
	write(10,*)"	real,parameter::Pforce = 1000.0"
	write(10,"(A27,E8.2)")"real,parameter::Emod = ", Emod
	write(10,"(A25,F4.1)")"real,parameter::nu = ", nu
	write(10,*)"	real,parameter::HookMat(3,3) = reshape(Emod/(1.0-nu**2)*(/ 1., nu, 0., nu, 1., 0., 0., 0., (1.0-nu)/2. /), (/3,3/))"
	write(10,*)
	write(10,*)"end module"
	close(10)
	
end subroutine