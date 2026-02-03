# Path of Compiler for My Macbook Air
#MPI_PATH = /usr/local
#INCLUDE = -I$(MPI_PATH)/include
#LIB = -L/usr/local/lib

#PRE_INCLUDE = -I/usr/local/metis/include
#PRE_LIB = -L/usr/local/lib -lmetis


# Path of Compiler for Dr. Ho's PC-Cluster in NCHC.
#MPI_PATH = /usr/mpi/gcc/openmpi-1.4.3
#INCLUDE = -I$(MPI_PATH)/include
#LIB = -L/usr/lib64

#PRE_INCLUDE = -I/home/c00nat00/CC_Folder/lib/include 
#PRE_LIB = -L/home/c00nat00/CC_Folder/lib/lib -lmetis


# Path of Compiler for ALPS PC-Cluster in NCHC.
MPI_PATH = /usr/mpi/intel/openmpi-1.8.4-qlc
INCLUDE = -I$(MPI_PATH)/include
LIB = -L/usr/lib64 -lmetis

PRE_INCLUDE = -I./ 
PRE_LIB = -L./ -lmetis


# Path of Compiler for Ariane.
#MPI_PATH = /opt/OpenMPI/gnu/1.6.2
#INCLUDE = -I$(MPI_PATH)/include
#LIB = -L/usr/lib64

#PRE_INCLUDE = -I/opt/metis/include
#PRE_LIB = -L/opt/metis/lib -lmetis



# C Compiler=====================================================
CC = g++ -O3 -c
CLINK = g++ -fPIC
#================================================================

# MPI Compiler===================================================
MPICC = $(MPI_PATH)/bin/mpicxx -O3 -c
LINK = $(MPI_PATH)/bin/mpicxx -fPIC
#================================================================

#================================================================
DEBUG = -g
#================================================================



# DSMC Solver
DSMC_FUNCTIONS_C = dsmc_main.o dsmc_class.o dsmc_init.o dsmc_function.o dsmc_output.o dsmc_readfile.o dsmc_toolfunction.o


# DSMC Preprocessor and Postprocessor
DSMC_PRE_FUNCTIONS = dsmc_main_preprocessor.o dsmc_class_pre.o dsmc_preprocessor.o
DSMC_POST_FUNCTIONS = dsmc_main_postprocessor.o dsmc_class.o dsmc_readfile.o dsmc_postprocessor.o


all : dsmc dsmc-pre dsmc-post

dsmc : ${DSMC_FUNCTIONS_C}
	${LINK} ${DEBUG} -o dsmc ${DSMC_FUNCTIONS_C}

dsmc-pre : ${DSMC_PRE_FUNCTIONS}
	${CLINK} ${DEBUG} -o dsmc-pre ${DSMC_PRE_FUNCTIONS} ${PRE_LIB}
	
dsmc-post : ${DSMC_POST_FUNCTIONS}
	${LINK} ${DEBUG} -o dsmc-post ${DSMC_POST_FUNCTIONS}

#================================================================

dsmc_main.o : dsmc_main.cpp
	${MPICC} ${DEBUG} dsmc_main.cpp ${INCLUDE}

dsmc_class.o : dsmc_class.cpp
	${MPICC} ${DEBUG} dsmc_class.cpp ${INCLUDE}

dsmc_init.o : dsmc_init.cpp
	${MPICC} ${DEBUG} dsmc_init.cpp ${INCLUDE}

dsmc_function.o : dsmc_function.cpp
	${MPICC} ${DEBUG} dsmc_function.cpp ${INCLUDE}

dsmc_output.o : dsmc_output.cpp
	${MPICC} ${DEBUG} dsmc_output.cpp ${INCLUDE}

dsmc_readfile.o : dsmc_readfile.cpp
	${MPICC} ${DEBUG} dsmc_readfile.cpp ${INCLUDE}

dsmc_toolfunction.o : dsmc_toolfunction.cpp 
	${MPICC} ${DEBUG} dsmc_toolfunction.cpp ${INCLUDE}

#================================================================

dsmc_main_preprocessor.o : dsmc_main_preprocessor.cpp
	${CC} ${DEBUG} dsmc_main_preprocessor.cpp ${PRE_INCLUDE}
	
dsmc_class_pre.o : dsmc_class_pre.cpp
	${CC} ${DEBUG} dsmc_class_pre.cpp ${PRE_INCLUDE}

dsmc_preprocessor.o : dsmc_preprocessor.cpp
	${CC} ${DEBUG} dsmc_preprocessor.cpp ${PRE_INCLUDE}

#================================================================

dsmc_main_postprocessor.o : dsmc_main_postprocessor.cpp
	${CC} ${DEBUG} dsmc_main_postprocessor.cpp ${INCLUDE}
	
dsmc_postprocessor.o: dsmc_postprocessor.cpp
	${CC} ${DEBUG} dsmc_postprocessor.cpp ${INCLUDE}

#================================================================

clean : 
	rm -f dsmc dsmc-pre dsmc-post *.o

clean_dat :
	rm -f *.dat

clean_inp :
	rm -f *.inp
