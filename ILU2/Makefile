############################################################
SYSHOME := $(shell /bin/pwd)
include $(SYSHOME)/Rules.make


############################################################
LIBOBJS = bcg.o iluoo.o ilu0.o dsort.o


############################################################
all:	    help

lib:	    $(LIBILU) info

clean:
	    @/bin/rm -f *.o *~ 

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib   - create library libilu.a"
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Libraries are located in lib/"
	    @echo " "


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o 


############################################################
$(LIBILU): $(LIBOBJS)
	   @echo ''
	   @rm -f $(LIBILU)
	   @ar -r $(LIBILU)  $(LIBOBJS)


############################################################

