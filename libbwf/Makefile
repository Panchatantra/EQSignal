CXX=g++
CXX_FLAGS=-I. -fPIC -shared -O3

SOURCE=$(wildcard *.cpp)
OBJECT=$(patsubst %.cpp, %.o, $(SOURCE))

ifeq ($(shell uname),Darwin)
	TARGET=libbwf.dylib
else
	TARGET=libbwf.so
endif

.SUFFIXES: .cpp .o


.cpp.o:
	$(CXX) -c $< $(CXX_FLAGS) -o $@

all: $(TARGET)

$(TARGET): $(OBJECT)
	$(CXX) $(OBJECT) $(CXX_FLAGS) -o $(TARGET)

install:$(TARGET)
	cp $(TARGET) ../EQSignal/lib
	cp $(TARGET) /usr/lib
	
clean:
	rm -vf $(OBJECT) $(TARGET)
