CXX=clang++
CXXFLAGS:=-Wall -Wextra -pedantic -std=c++20
RELEASEFLAGS:=-Ofast

OUTPUT := trains
BONUS := trains-bonus

.PHONY: all clean bonus

all: $(OUTPUT) 

$(OUTPUT): simulate.cc main.cc
	mpicxx $(CXXFLAGS) $(RELEASEFLAGS) -o $@ $^
	
clean:
	$(RM) *.o $(OUTPUT)

bonus: $(BONUS)

$(BONUS): simulate_bonus.cc main_bonus.cc
	mpicxx $(CXXFLAGS) $(RELEASEFLAGS) -o $@ $^
