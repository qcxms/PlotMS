 PROG = ./plotms 
#--------------------------------------------------------------------------
 OSTYPE=LINUXI
#--------------------------------------------------------------------------

OBJS1=plotms_v5.0.o

OBJS2 = 

OBJS = $(OBJS1) $(OBJS2)
#--------------------------------------------------------------------------

ifeq ($(OSTYPE),LINUXI)
  PREOPTS =
  FC = ifort
  CC = gcc
  LINKER = ifort
  LIBS    = 
  PREFLAG = -E -P
  CCFLAGS = -O -DLINUX 
  FFLAGS = -O2 
#-w90 -tpp6 -d0 -C
endif                     

ifeq ($(OSTYPE),LINUXL)
  FC = lf95
  CC = gcc
  LINKER = lf95 --staticlink
  PREFLAG = -E -P
  CCFLAGS = -O -DLINUX
# FFLAGS = -O --ntrace --tpp --info 
  FFLAGS = -O --trace --tpp --info --chk a,e,s,u --prefetch 2
endif


# diese ziele gibts:
.PHONY: all
.PHONY: clean
# dieses ist das erste auftretende,
# wird also beim aufruf von make erzeugt (default)
all: $(PROG)


#--------------------------------------------------------------------------
# example.f: printversion.h
# example.f haengt von printversion.h ab
# wenn sich also  printversion.h aendert, wird example.f
# (und damit auch example.o) neu gemacht.
# was auch geht:

#--------------------------------------------------------------------------
# implizite Regel zur Erzeugung von *.o aus *.F ausschalten
%.o: %.F

# aus *.F mache ein *.f
%.f: %.F
	@echo "making $@ from $<"
	$(CC) $(PREFLAG) $(PREOPTS) $< -o $@

# aus *.f mache ein *.o
%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# aus *.c mache ein *.o
%.o: %.c
	@echo "making $@ from $<"
	$(CC) $(CCFLAGS) -c $< -o $@

# linken
$(PROG): $(OBJS) 
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)


#aufraeumen
clean:
	rm -f *.o $(PROG) 
	rm -f $(patsubst %.F, %.f, $(wildcard *.F))


