
name := udod

# External libraries:
# ATTENTION: do not use cuba-3.0! it does not work for us
#------------------------------------------------------------------------------
# CLHEP
#------------------------------------------------------------------------------
#CLHEP_DIR := /opt/geant4/clhep/2.0.4.2/
#CLHEPINC  := -I${CLHEP_DIR}/include
#CLHEPLIB  := ${CLHEP_DIR}/lib/libCLHEP.a
#------------------------------------------------------------------------------
# HEPMC
#------------------------------------------------------------------------------
# HEPMC_DIR := /home/afs/cern.ch/sw/lcg/external/HepMC/2.06.09/i686-deb-gcc47
#HEPMC_DIR := /opt/HepMC-2.06.04/
#HEPMCINC  := -I${HEPMC_DIR}/include
#HEPMCLIB  := ${HEPMC_DIR}/lib/libHEPMC.a
#------------------------------------------------------------------------------
# ROOT
#------------------------------------------------------------------------------
ROOTCFLAGS := $(shell root-config --cflags) 
#ROOTLIBS   := $(shell root-config --libs)
#ROOTGLIBS  := $(shell root-config --glibs)
ROOTLIBS = -L/opt/root/6.10.08/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic

#------------------------------------------------------------------------------

# Compiler flags
#------------------------------------------------------------------------------
CPP        := g++
CPPFLAGS   := -g -Wall -pipe -fPIC
CPPFLAGS   += -MMD  # generate dependency files (not system headers)

# includes for libs 
CPPFLAGS   += -I./Udod -I./cuba_int -I./Cuba-2.1 -I./models
CPPFLAGS   += ${CLHEPINC}
CPPFLAGS   += ${HEPMCINC}
CPPFLAGS   += ${ROOTCFLAGS}

# do not use -l cuba to prevent of using the system cuba-3 library
LDFLAGS    := -L./ -lUdod -lcuba_int ./libCuba.a 
#LDFLAGS    += -L${CLHEP_DIR}/lib -lCLHEP
#LDFLAGS    += -L${HEPMC_DIR}/lib -lHepMC
LDFLAGS    += ${ROOTLIBS}

#------------------------------------------------------------------------------
object_dir      := $(shell uname)
PROGRAM         := $(name).exe
LIBUDOD         := libUdod.a
LIBCUBA_INT     := libcuba_int.a

vegas_srcs      := $(wildcard Vegas/*.cxx)
vegas_objs      := $(patsubst Vegas/%.cxx, %.o, $(vegas_srcs))
vegas_objs      := $(addprefix $(object_dir)/, $(vegas_objs))

udod_srcs       := $(wildcard  Udod/*.cxx)
udod_objs       := $(patsubst Udod/%.cxx, %.o, $(udod_srcs))
udod_objs       := $(addprefix $(object_dir)/, $(udod_objs))

cuba_int_srcs   := $(wildcard cuba_int/*.cxx)
cuba_int_objs   := $(patsubst cuba_int/%.cxx, %.o, $(cuba_int_srcs))
cuba_int_objs   := $(addprefix $(object_dir)/, ${cuba_int_objs})

models_srcs       := $(wildcard  models/*.cxx)
models_objs       := $(patsubst models/%.cxx, %.o, $(models_srcs))
models_objs       := $(addprefix $(object_dir)/, $(models_objs))

VPATH           := cuba_int/ Udod/ Vegas/ models/
SUFFIXES        := .o .cxx
#------------------------------------------------------------------------------
# rules
$(object_dir)/%.o: %.cxx
	${CPP} ${CPPFLAGS} -c -o $@ $<
#------------------------------------------------------------------------------

.PHONY: all clean distclean memcheck

all: $(object_dir) $(PROGRAM)

$(PROGRAM): libCuba.a $(LIBCUBA_INT) $(LIBUDOD)
	${CPP} main.cxx ${CPPFLAGS} ${LDFLAGS} -o $(PROGRAM)

libCuba.a:
	$(MAKE) --makefile=Cuba-2.1/Makefile srcdir=Cuba-2.1

$(LIBUDOD): $(udod_objs) $(models_objs)
	$(AR) $(ARFLAGS) $@ $?
	ranlib $@

$(LIBCUBA_INT): $(cuba_int_objs)
	$(AR) $(ARFLAGS) $@ $?
	ranlib $@

libVegas.a: $(vegas_objs)
	$(AR) $(ARFLAGS) $@ $?
	ranlib $@

libUdod.so: $(udod_objs) $(models_objs)
	${CPP} -shared $^ -o$@

libcuba_int.so: $(cuba_int_objs)
	${CPP} -shared $^ -o$@

$(object_dir):
	@mkdir -p $(object_dir)

clean:
	rm -fr ./$(object_dir) *.d test_cubaint test_cktable
	rm -fr core
	find . -name "*~" -exec rm -fr {} \;

distclean: clean
	rm -f lib*.a *.so $(PROGRAM) *.root

memcheck: $(PROGRAM)
	valgrind --workaround-gcc296-bugs=yes -v --leak-check=full ./$(PROGRAM)

# Debugging Makefiles: http://www.ddj.com/development-tools/197003338
print-%: ; @echo $* is $($*)

# testing
test_cubaint: 
	${CPP} cuba_int/cuba_int.cxx -DMAINCUBAINT ${CPPFLAGS} \
                                ${CUBALIB} -o $@

test_cktable: 
	${CPP} Udod/CkTable.cxx -DMAINCKTABLE ${CPPFLAGS} -o $@

# dependencies (.d-files are created by gcc at compile time)
INCDEP := $(wildcard *.d)
INCDEP += $(wildcard $(object_dir)/*.d)
ifneq ($(INCDEP),)
include $(INCDEP)
endif

