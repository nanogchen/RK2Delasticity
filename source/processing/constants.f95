 module constants
 	implicit none

 	!domain geometry
   	real,parameter::x_start =  0.0
   	real,parameter::x_end =  4.0
   	real,parameter::y_start = -0.5
   	real,parameter::y_end =  0.5
 	real,parameter::Length = x_end - x_start
 	real,parameter::Height = y_end - y_start

 	!domain node/Guass point generating parameter
    integer,parameter::npx =  41
    integer,parameter::npy =  11
 	integer,parameter::np = npx*npy
    integer,parameter::nGx =   2
    integer,parameter::nGy =   2
 	integer,parameter::npxGauss = nGx*(npx-1)
 	integer,parameter::npyGauss = nGy*(npy-1)
 	integer(kind=4),parameter::npGauss = npxGauss*npyGauss

 	!boundary Guass point generating parameter
    integer,parameter::nGBC =   2
 	integer,parameter::npBC = 2*(npx+npy)-4
 	integer,parameter::npInner = np - npBC
    integer,parameter::npEBC =  11
    integer,parameter::npNBC =  11
    integer,parameter::npGaussEBC =  20
    integer,parameter::npGaussNBC =  20

 	!parameter for penalty method
    real,parameter::beta = 0.10E+06
    integer,parameter::i_algorithm =   1

 	!Nodal configuration
    real,parameter::dilation = 1.75
 	real,parameter::dx = (x_end - x_start)/(npx-1)
 	real,parameter::dy = (y_end - y_start)/(npy-1)
 	real,parameter::hx = dilation*dx ! for rectangular phi
 	real,parameter::hy = dilation*dy ! for rectangular phi
 	real,parameter::hr = dilation*sqrt(dy**2+dx**2) ! for circular phi

 	!Neighbor list to search neighbor Gauss points for each nodes
 	integer,parameter::NumNeighborInner = (( int(dilation) + 1 )*2)**2*nGx*nGy
 	integer,parameter::NumNeighborEBC = 16*nGBC ! the maximum num of neighbors for an essential boundary node
 	integer,parameter::NumNeighborNBC = 16*nGBC ! the maximum num of neighbors for a natural boundary node

 	!material properties
 	real,parameter::Imoment = Height**3/12.0
 	real,parameter::Pforce = 0.10E+04
    real,parameter::Emod = 0.30E+08
    real,parameter::nu =  0.3
   integer,parameter::i_model =   1
 	real,parameter::HookMat_strain(3,3) = reshape(Emod*(1-nu)/((1+nu)*(1-2*nu))*(/ 1., nu/(1-nu), 0., & 
 				nu/(1-nu), 1., 0., 0., 0., (1.0-2*nu)/(2.*(1-nu)) /), (/3,3/))
 	real,parameter::HookMat_stress(3,3) = reshape(Emod/(1.0-nu**2)*(/ 1., nu, 0., nu, 1., 0., 0., 0., (1.0-nu)/2. /), (/3,3/))

 end module
