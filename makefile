CXX = g++

TARGET = main

SRC = main.cpp Table.cpp Simulation.cpp
OBJ = $(SRC:.cpp=.o)
DEP = $(OBJ:.o=.d)

PGO_DIR = ./pgo

BASE_CXXFLAGS = -Ofast -march=native -mtune=native -fno-math-errno \
                -std=c++23 -fopenmp -flto=auto \
                -Wall -I./include -MMD -MP

BASE_LDFLAGS = -Ofast -march=native -mtune=native \
               -fopenmp -flto=auto

CXXFLAGS = $(BASE_CXXFLAGS)
LDFLAGS  = $(BASE_LDFLAGS)

ifeq ($(PGO),gen)
CXXFLAGS += -fprofile-generate=$(PGO_DIR) -fprofile-update=atomic
LDFLAGS  += -fprofile-generate=$(PGO_DIR) -fprofile-update=atomic
endif

ifeq ($(PGO),use)
CXXFLAGS += -fprofile-use=$(PGO_DIR) -fprofile-correction
LDFLAGS  += -fprofile-use=$(PGO_DIR) -fprofile-correction
endif

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include $(DEP)

pgo-gen:
	rm -rf $(PGO_DIR)
	mkdir -p $(PGO_DIR)
	$(MAKE) clean-build
	$(MAKE) PGO=gen all

pgo-use:
	$(MAKE) clean-build
	$(MAKE) PGO=use all

clean-build:
	rm -f $(OBJ) $(DEP) $(TARGET)

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)
	rm -rf $(PGO_DIR)

rebuild: clean all

.PHONY: all clean clean-build rebuild pgo-gen pgo-use