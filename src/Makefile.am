INCLUDES = -I/opt/local/include
bin_PROGRAMS = SlaterGen
SlaterGen_SOURCES = main.cc mpreal.cpp dlmalloc.c root_weight.cc util.cc generate.cc
SlaterGen_LDADD = $(FLIBS) -L/opt/local/lib -lmpfr -lgmp -lblas
ACLOCAL_AMFLAGS = -I m4
