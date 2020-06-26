IDIR=${CURDIR}/src

	
program=novaSTA

SOURCES := $(shell find $(IDIR) -name '*.cpp')
SOURCES +=$(shell find $(IDIR) -name '*.cu')
OBJECTS = $(SOURCES:.cpp=.o)
OBJECTS += $(SOURCES:.cu=.o) 

CXX = mpic++ -O3 -s -DNDEBUG
CXXFLAGS = -std=c++14 #-std=c++0x

CUDA = nvcc
CUDA_INCLUDES = -I"${cudapath}/include"

#OBJECTS = $(SOURCES:.cpp=.o)

LIBRARIES = -lfftw3 -lfftw3f -lpthread -lcudart -lcufft

all: ${program}

build: ${program}

debug: CXX = mpic++ -g -W #-Wall -Werror 
debug: LDFLAGS += -g
debug: ${program}

 
${program}: CXXFLAGS += $(foreach d, $(includepath), -I$d)
${program}: LDFLAGS += $(foreach d, $(libpath), -L$d)

${program}: $(OBJECTS)
	$(CXX) -o $@ $(LDFLAGS) $(IDIR)/*.o $(LIBRARIES)

clean:
	rm -f src/*.o ${program} src/*.d

.cpp.o:
	$(CXX) -MD -MP $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)

%.o: %.cu 
	$(CUDA) -c $(CUDA_INCLUDES) $< -o $@

