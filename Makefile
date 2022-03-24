# 1. Put this file in the same folder as your 'driver' code 
#    (the code containing the 'main' function).

# 2. Edit LIBRARY_DIR to point at the location of your ITensor Library
#    source folder (this is the folder that has options.mk in it)
LIBRARY_DIR=/home/shong/ITensor-ITensor-ac434f2/

# 3. If your 'main' function is in a file called 'myappname.cc', then
#    set APP to 'myappname'. Running 'make' will compile the app.
#    Running 'make debug' will make a program called 'myappname-g'
#    which includes debugging symbols and can be used in gdb (Gnu debugger);
APP=TRG_it

# 4. Add any headers your program depends on here. The make program
#    will auto-detect if these headers have changed and recompile your app.
HEADERS=complex_bessel.h complex_def.h


# 5. For any additional .cc files making up your project,
#    add their full filenames here.
CCFILES=$(APP).cc 

#################################################################
#################################################################
#################################################################
#################################################################


include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/all.h

#Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $< -L/home/shong/apps/complex_bessel-release-0.6/build -lgfortran -lcomplex_bessel

.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $< -L/home/shong/apps/complex_bessel-release-0.6/build -lgfortran -lcomplex_bessel


#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS) -L/home/shong/apps/complex_bessel-release-0.6/build -lgfortran -lcomplex_bessel

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS) -L/home/shong/apps/complex_bessel-release-0.6/build -lgfortran -lcomplex_bessel


clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g

mkdebugdir:
	mkdir -p .debug_objs

