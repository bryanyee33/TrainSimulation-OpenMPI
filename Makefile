CXX=clang++
CXXFLAGS:=-Wall -Wextra -pedantic -std=c++20
RELEASEFLAGS:=-Ofast

OUTPUT := trains

.PHONY: all clean

all: $(OUTPUT) 

$(OUTPUT): simulate.cc main.cc
	mpicxx $(CXXFLAGS) $(RELEASEFLAGS) -o $@ $^
	
clean:
	$(RM) *.o $(OUTPUT)
