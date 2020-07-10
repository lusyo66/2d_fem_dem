EIGEN_ROOT = /usr/local/eigen-3.3.4
MPICXX = mpic++
OPENMP = -fopenmp
C++11 - -std=c++11

CXXFLAGS = $(OPENMP)
EIGEN_INCLUDE=-I$(EIGEN_ROOT)
INCLUDE=$(EIGEN_INCLUDE)

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

EXECUTABLE = 2dFEM

$(EXECUTABLE): $(OBJECTS)
	$(MPICXX) -o $@ $(CXXFLAGS) $(OBJECTS)

%.o: %.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDE)

sinclude $(SOURCES:.cpp=.d)

%.d: %.cpp
	$(MPICXX) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

clean:
	-rm -f *.o *.d  *~ *.tar $(EXECUTABLE)


