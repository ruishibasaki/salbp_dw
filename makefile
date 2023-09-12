# compiler
CC = g++

# module name
NAME = salbp

# basic directory
DIR = ./

# debug switches
#SW = -Wall -ggdb3
# production switches
SW = -w -O3

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(NAME)


# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)*.o $(DIR)*~ $(NAME)

# define include the necessary libs- - - - - - - - - - - - - - - - - - -

SYSTEM = x86-64_linux
LIBFORMAT = static_pic
CONCERTDIR = /opt/ibm/ILOG/CPLEX_Studio1210/concert
CPLEXDIR = /opt/ibm/ILOG/CPLEX_Studio1210/cplex

#SYSTEM = x86-64_linux
#LIBFORMAT = static_pic
#CONCERTDIR = /users/rsashibasaki/ibm/ILOG/CPLEX_Studio1210/concert
#CPLEXDIR = /users/rsashibasaki/ibm/ILOG/CPLEX_Studio1210/cplex

#SYSTEM = x86-64_osx
#LIBFORMAT = static_pic
#CONCERTDIR = /Applications/CPLEX_Studio128/concert
#CPLEXDIR = /Applications/CPLEX_Studio128/cplex


CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXINC = -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include


# ---------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------
CPLEXFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -ldl -lpthread 




OBJ = pattern_manager.o solver.o main.o data.o heuristic.o

./heuristic.o: heuristic.cpp
	$(CC) -c $< -o $@ $(CPLEXINC) $(SW)
./solver.o: solver.cpp
	$(CC) -c $< -o $@ $(CPLEXINC) $(SW) 
./main.o: main.cpp
	$(CC) -c $< -o $@ $(CPLEXINC) $(SW) 
./pattern_manager.o: pattern_manager.cpp
	$(CC) -c $< -o $@ $(SW) 
./data.o: data.cpp
	$(CC) -c $< -o $@ $(SW) 


$(NAME): $(OBJ)
	$(CC) -o  $(NAME) $(OBJ) $(CPLEXFLAGS)  $(SW) 
