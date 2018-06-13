#------------------------------------------------------------------------------
# CFLAGS  := -g -Wall -pipe
CFLAGS  := -O2 -pipe -fPIC
# CFLAGS  += -MMD  # generate dependency files (not system headers)

ARFLAGS := rU

RANLIB  := ranlib
#------------------------------------------------------------------------------
srcdir  := .
LIB     := libCuba.a

common = $(srcdir)/src/common
vegas = $(srcdir)/src/vegas
suave = $(srcdir)/src/suave
divonne = $(srcdir)/src/divonne
cuhre = $(srcdir)/src/cuhre

CFLAGS  += -I$(srcdir) -I$(common)

COMMON_DEPS = $(common)/stddecl.h $(common)/ChiSquare.c $(common)/Erf.c

VEGAS_C = $(vegas)/Vegas.c
VEGAS_DEPS = $(COMMON_DEPS) $(common)/Random.c \
  $(vegas)/decl.h $(vegas)/common.c \
  $(vegas)/Grid.c $(vegas)/Integrate.c

SUAVE_C = $(suave)/Suave.c
SUAVE_DEPS = $(COMMON_DEPS) $(common)/Random.c \
  $(suave)/decl.h $(suave)/common.c \
  $(suave)/Grid.c $(suave)/Fluct.c $(suave)/Sample.c \
  $(suave)/Integrate.c

DIVONNE_C = $(divonne)/Divonne.c
DIVONNE_DEPS = $(COMMON_DEPS) $(common)/Random.c \
  $(divonne)/decl.h $(divonne)/common.c \
  $(divonne)/KorobovCoeff.c $(divonne)/Rule.c $(divonne)/Sample.c \
  $(divonne)/FindMinimum.c $(divonne)/Explore.c $(divonne)/Split.c \
  $(divonne)/Integrate.c

CUHRE_C = $(cuhre)/Cuhre.c
CUHRE_DEPS = $(COMMON_DEPS) \
  $(cuhre)/decl.h $(cuhre)/common.c \
  $(cuhre)/Rule.c $(cuhre)/Integrate.c
#------------------------------------------------------------------------------

.PHONY: all clean

all: $(LIB)

$(LIB): $(LIB)(Vegas.o) $(LIB)(Suave.o) $(LIB)(Divonne.o) $(LIB)(Cuhre.o)
	-$(RANLIB) $(LIB)

$(LIB)(Vegas.o): $(VEGAS_C) $(VEGAS_DEPS)
	$(CC) $(CFLAGS) -I$(vegas) -c -o Vegas.o $(VEGAS_C)
	$(AR) $(ARFLAGS) $(LIB) Vegas.o
	$(RM) Vegas.o

$(LIB)(Suave.o): $(SUAVE_C) $(SUAVE_DEPS)
	$(CC) $(CFLAGS) -I$(suave) -c -o Suave.o $(SUAVE_C)
	$(AR) $(ARFLAGS) $(LIB) Suave.o
	$(RM) Suave.o

$(LIB)(Divonne.o): $(DIVONNE_C) $(DIVONNE_DEPS)
	$(CC) $(CFLAGS) -I$(divonne) -c -o Divonne.o $(DIVONNE_C)
	$(AR) $(ARFLAGS) $(LIB) Divonne.o
	$(RM) Divonne.o

$(LIB)(Cuhre.o): $(CUHRE_C) $(CUHRE_DEPS)
	$(CC) $(CFLAGS) -I$(cuhre) -c -o Cuhre.o $(CUHRE_C)
	$(AR) $(ARFLAGS) $(LIB) Cuhre.o
	$(RM) Cuhre.o

clean:
	$(RM) $(LIB)

# Debugging Makefiles: http://www.ddj.com/development-tools/197003338
print-%: ; @echo $* is $($*)

# dependencies (.d-files are created by gcc at compile time)
# INCDEP := $(wildcard *.d)
# INCDEP += $(wildcard $(object_dir)/*.d)
# ifneq ($(INCDEP),)
# include $(INCDEP)
# endif

