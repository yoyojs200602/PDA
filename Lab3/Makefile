main_obj = main.cpp
INC = abacus/abacus.hpp FM/FM.hpp

GXX = g++
FLAG = -std=c++17 -pthread -g -w -O3
EXE = Lab3

Lab3: main.cpp $(INC)
	$(GXX) $(FLAG) -I $(INC) -o $(EXE) $(main_obj)

clean:
	rm Lab3
