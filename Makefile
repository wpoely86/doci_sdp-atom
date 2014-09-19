###############################################################################
#
#  makefile template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Sources for all modules
# -----------------------------------------------------------------------------
BINNAME = doci_sdp
CPPSRC	= doci_sdp.cpp\
            Matrix.cpp\
	    Vector.cpp\
	    BlockStructure.cpp\
	    Container.cpp\
	    helpers.cpp\
	    Tools.cpp\
	    TPM.cpp\
	    SPM.cpp\
	    SUP.cpp\
	    EIG.cpp\
	    Lineq.cpp\
#            SPM.cpp\
#            TPM.cpp\
#            PHM.cpp\
#            DPM.cpp\
#            PPHM.cpp\
#            SUP.cpp\
#            EIG.cpp\

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

BRIGHT_ROOT= .

INCLUDE = ./include  -I../integrals/include -I/usr/include/libint2

LIBS= -llapack -lblas -lhdf5 integrals.a -lgsl -lint2

CC	= gcc
CXX	= clang++

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS	= -I$(INCLUDE) -std=c++11 -g -Wall -O2 -march=native -Wno-unused-variable
LDFLAGS	= -g -Wall -O2 


# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	@echo
	@echo '  +++ Building $(BINNAME)...'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQG"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built successfully!'; \
	   echo; \
	 fi

#------------------------------------------------------------------------------
#  Compile with only P and Q conditions activated
#------------------------------------------------------------------------------

PQ:
	@echo
	@echo '  +++ Building $(BINNAME) with P and Q conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQ"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P and Q conditions successfully!'; \
	   echo; \
	 fi

PQG:
	@echo
	@echo '  +++ Building $(BINNAME) with P, Q and G conditions active'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQG"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P, Q and G conditions successfully!'; \
	   echo; \
	 fi

PQGT1:
	@echo
	@echo '  +++ Building $(BINNAME) with P, Q ,G and T1 conditions active'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQGT1"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P, Q, G and T1 conditions successfully!'; \
	   echo; \
	 fi

PQGT:
	@echo
	@echo '  +++ Building $(BINNAME) with P, Q ,G and T1 conditions active'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQGT"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P, Q, G, T1 and T2 conditions successfully!'; \
	   echo; \
	 fi

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for Makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) -c $(FFLAGS) $(SFLAGS) $(@:.o=.for) -o $@

%.o:	%.c Makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) -c $(CFLAGS) $(SFLAGS) $(@:.o=.c) -o $@

%.o:	%.cpp Makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) -c $(CFLAGS) $(SFLAGS) $(DEFS) $(@:.o=.cpp) -o $@


# -----------------------------------------------------------------------------
#   Link everything together
# -----------------------------------------------------------------------------
$(BRIGHT_ROOT)/$(BINNAME):	Makefile $(OBJ) 
	@echo; echo "Linker: creating $(BRIGHT_ROOT)/$(BINNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/$(BINNAME) $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all object files
# -----------------------------------------------------------------------------
clean:
	@echo -n '  +++ Cleaning all object files ... '
	@echo -n $(OBJ)
	@rm -f $(OBJ)
	@echo 'Done.'

# -----------------------------------------------------------------------------
#   Make new documentation using doxygen
# -----------------------------------------------------------------------------
doc:
	@doxygen doc-config

# ====================== End of file 'makefile.in' ========================== #
