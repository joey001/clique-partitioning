
#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

TARGET = cpp
DIR_SRC = ./src
DIR_OBJ = ./
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0
# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O3 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

SRC = $(wildcard ${DIR_SRC}/*.cpp) 
OBJ = $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC})) 
#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
all:
	make $(TARGET)
	/bin/rm $(OBJ)


# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ 
	/bin/rm -rf $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ALERT : to link the objectives, (OBJ) must put before CCLNDIRS and CCLNFLAGS
# ALERT : The indent is must be a TAB (\t), not 4 blanks
${TARGET}: ${OBJ}
	$(CCC) $(OBJ) -o $@


${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
	$(CCC) -c $(CCOPT) $< -o $@

# Local Variables:
# mode: makefile
# End:
