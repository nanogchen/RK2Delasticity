FC = gfortran

prog = RKpost.exe

# module file end with .f95 and .f03
modsrcs1 = $(wildcard *.f95)

srcs = $(wildcard *.f90)
		
objs = $(patsubst %.f90, %.o, $(srcs)) $(patsubst %.f95, %.o, $(modsrcs1))

mods = $(patsubst %.f95, %.mod, $(modsrcs1))

$(prog): $(objs)
	$(FC) -g -fcheck=all -Wall -o $(prog) $(objs)	
	
$(objs): $(srcs) $(mods)
	$(FC) -g -fcheck=all -c $(srcs)	
	
$(mods): $(modsrcs1) 
	$(FC) -g -fcheck=all -c $(modsrcs1)
	
clean:
	rm $(prog) $(objs) $(mods)