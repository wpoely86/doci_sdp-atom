###############################################################################
#
#  makefile template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Sources for all modules
# -----------------------------------------------------------------------------
# prefix for binaries
BINNAME = doci

CPPSRC	=   Matrix.cpp\
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
            PHM.cpp\
	    BoundaryPoint.cpp\
	    PotentialReduction.cpp\
	    SimulatedAnnealing.cpp\
	    LocalMinimizer.cpp\
#            DPM.cpp\
#            PPHM.cpp\

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

BRIGHT_ROOT= .

INCLUDE = -Iinclude -Iextern/include

LIBS= -llapack -lblas -lhdf5 -Lextern -lsimanneal
# for the MKL:
#LIBS= -lhdf5 -Lextern -lsimanneal -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

ifeq ($(origin CC), default)
    CC = clang
endif

ifeq ($(origin CXX), default)
    CXX = mpicxx
endif

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS	= $(INCLUDE) -std=c++11 -g -Wall -O2 -march=native -Wno-unknown-pragmas -Wno-sign-compare
LDFLAGS	= -g -Wall -O2


# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	$(MAKE) -C extern
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

P:
	@echo
	@echo '  +++ Building $(BINNAME) with P condition'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS=""
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P condition successfully!'; \
	   echo; \
	 fi

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
$(BRIGHT_ROOT)/$(BINNAME):	Makefile $(OBJ) doci.o
	@echo; echo "Linker: creating $(BRIGHT_ROOT)/$(BINNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/$(BINNAME) doci.o $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all object files
# -----------------------------------------------------------------------------
clean:
	$(MAKE) -C extern clean
	@echo -n '  +++ Cleaning all object files ... '
	@echo -n $(OBJ)
	@rm -f $(OBJ) doci_bp.o doci_sdp.o
	@echo 'Done.'

# -----------------------------------------------------------------------------
#   Make new documentation using doxygen
# -----------------------------------------------------------------------------
doc:
	@doxygen doc-config

bp: $(OBJ) doci_bp.o
	@echo 'Building boundary point method'
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/doci_bp doci_bp.o $(OBJ) $(LIBS)

sdp: $(OBJ) doci_sdp.o
	@echo 'Building Potential reduction method'
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/doci_sdp doci_sdp.o $(OBJ) $(LIBS)

print: $(OBJ) print.o
	@echo 'Building print'
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/print print.o $(OBJ) $(LIBS)

# ====================== End of file 'makefile.in' ========================== #
