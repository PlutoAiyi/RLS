# Makefile

CXX = g++

CXXFLAGS = -O3 -g

TARGET = RLS 

SRC = main.cpp Graph.cpp Utility.cpp RandList.hpp

DEPS = Graph.h Utility.h

all: $(TARGET)

$(TARGET): $(SRC) $(DEPS)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all clean
