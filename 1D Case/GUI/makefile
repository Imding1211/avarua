
all: $(OBJECTS) 
	f2py -c ATM_2case_26method_solver.f -m solver

GUI:
	python ATM_2case_26method_main_tk.py
	
.PHONY:clean
clear:  
	rm -f *.exe
	rm -f *.mod
	rm -f *.o
	rm -f *.dat
	rm -f *.out