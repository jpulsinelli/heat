MACHINE = summit_xlf

OPT_LEVEL = OPTIMIZE
FLAGS     = $(FLAGS_$(OPT_LEVEL))

include ./Makefile_Main

all: Heat

Heat: \
	$(objects) Heat.o
	$(FLINKER) $(FLAGS) -o Heat_$(MACHINE) \
	$(objects) Heat.o $(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Heat_$(MACHINE)
