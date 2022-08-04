
all: $(OBJECTS) 
	f2py -c ATM_2case_6method_solver.f -m solver
	python ATM_6method_main.py

ref:
	f2py -c ATM_2case_6method_solver.f -m solver
	python 12800.py

run:
	f2py -c ATM_2case_6method_solver.f -m solver
	nohup python ATM_6method_main.py &
	
.PHONY:clean
clear:  
	rm -f *.exe
	rm -f *.mod
	rm -f *.o
	rm -f *.dat
	rm -f *.out