
all: $(OBJECTS) 
	f2py -c PDE_solver.f90 -m solver
	python PDE_main.py

run:
	f2py -c PDE_solver.f90 -m solver
	nohup python PDE_main.py &
	
.PHONY:clean
clear:  
	rm -f *.exe
	rm -f *.mod
	rm -f *.o
	rm -f *.dat
	rm -f *.out