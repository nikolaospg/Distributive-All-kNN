# make (all) command: 	compiles and produces executables for all of the implementations
# make run command: 	runs the executable of the MPI implementations, with the proper inputs
# make clean command:	removes all the executable files created by the make command

CC=gcc
MPICC=mpicc

default: all

V0: V0.c
	$(CC) -o V0 V0.c -lopenblas

V1: V1.c
	$(MPICC) -o V1 V1.c -lopenblas -lm

V2:	V2.c
	$(MPICC) -o V2 V2.c -lopenblas -lm

all: V0 V1 V2


# If the first argument is "run"...
ifeq (run, $(firstword $(MAKECMDGOALS)))

  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2, $(words $(MAKECMDGOALS)), $(MAKECMDGOALS))

  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif


### There are 2 ways for the user to use make run: ###

# 1) The user can select the file to be read and the variable k.
# e.g. make run processes=2 [path/to/matrix/features.csv] 50
# That means that the ./V0, ./V1 and ./V2 programs will read the matrix features.csv and run it with 2 processes for 50 nearest neighbours.
#
# 2) The user can select to run a random array with specific n, d and k.
# e.g. make run processes=2 1000 20 50
# That means that the ./V0, ./V1 and ./V2 programs will make a random matrix with n = 1000, d = 20 and k = 50 with 2 processes.
run:
	./V0 $(RUN_ARGS)
	mpirun -np $(processes) ./V1 $(RUN_ARGS)
	mpirun -np $(processes) ./V2 $(RUN_ARGS)

.PHONY: clean

clean:
	rm -f V0 V1 V2
