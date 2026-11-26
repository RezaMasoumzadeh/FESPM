.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

include Makefile.build

OBJECTS = \
geo_in3mg.o \
geo_load3mg.o \
geo_mark3mg.o \
geo_move3mg.o \
geo_gaus3mg.o \
geo_heat3mg.o \
geo_surf3mg.o \
geo_i3vis_output_paraview.o \
spmgrid.o \
spmupdate.o \
fespm.o \
spmsl.o


OBJECTS2 = \
geo_i3mg.o \
geo_load3mg.o \
geo_mark3mg.o \
geo_move3mg.o \
geo_gaus3mg.o \
geo_heat3mg.o \
geo_surf3mg.o \
geo_i3vis_output_paraview.o \
spmgrid.o \
spmupdate.o \
fespm.o \
spmsl.o

# These are always built
.PHONY: all clean resclean

# Default. Builds everything if outdated
all: in3mg i3mg 


.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

.f.o:
	$(F90) $(FLAGS) $*.f

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) $*.c



in3mg:	$(OBJECTS)
	$(CC) $(OPTIONS) $(OBJECTS)  $(LIBS)  -o in3mg 

i3mg:	$(OBJECTS2)
	$(CC) $(OPTIONS) $(OBJECTS2)  $(LIBS) -o i3mg 
	

clean:
	-rm -f  *.o   i3mg in3mg

	
# cleans results only	
resclean:
	-rm -f *.pvd *.prn *.vtr *.vtk *.cpf times.txt


