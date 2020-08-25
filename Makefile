OBJ=  ratios.o wiser.o f2nmc_new.o 

# -lstdc++ for Gagik's stuff
#
FOR   =  -lstdc++ -lnsl -lcrypt -ldl -lgfortran
CERNLIBS =  -L/apps/cernlib/x86_64_rhel6_4.7.2/2005/lib -lmathlib  -lpacklib
claspyth : $(OBJ)
	 gfortran  -o	ratios  $(OBJ) $(CERNLIBS) $(FOR) 
$(OBJ) : %.o: %.f
	gfortran   -DLinux -fno-automatic  -fno-second-underscore  -I. -I./include -c $< -o $@  

clean:
	rm -f ratios   $(OBJ)



##g77  -O2 -fno-automatic -finit-local-zero -ffixed-line-length-none -fno-second-underscore \
##        -DLinux \
##        -I. -I./ -I/group/clas/builds/release-4-14/packages/include -I/group/clas/builds/release-4-14/packages/inc_derived -I/apps/tcl/include -I/usr/X11R6/include -c \
##        aac.F  -o /home/avakian/w6/tmp/obj/LinuxRHEL3/generator/aac.o


















