# SWIG
INTERFACES = hello.i  orbCalc.i orbCalc2.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cxx)
PROXIES    = $(INTERFACES:.i=.py      )
# Default target: build the tests

.PHONY : all
all: cwrapp  
cwrapp: $(WRAPPERS) hello.cc orbCalc.cc orbCalc2.cc
	./setup.py build_ext -i

# Rule: %.i -> %_wrap.cpp
%_wrap.cxx: %.i %.h ./numpy.i
	swig -c++ -python $<


# Clean target
.PHONY : clean
clean:
	$(RM) -r build
	$(RM) *.so *.pyc *_wrap.h
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
	$(RM) *.o
