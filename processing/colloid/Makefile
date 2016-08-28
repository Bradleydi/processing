CC = g++

CXXFLAGS=-O2 -fPIC -Wall -pedantic -c 

MODULES=colloid_base \
		data_preprocess \
		dct \
		dynamics \
		elasticity \
		ellipsoid_normalmode \
		fftw \
		io \
		lapack \
		miscellaneous \
		normalmode \
		plot \
		qchi4  \
		statistics \
		structure 
SRC_DIR=$(addprefix src/,$(MODULES))
BUILD_DIR=$(addprefix build/,$(MODULES))

SRC=$(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
OBJ= $(patsubst src/%.cpp,build/%.o,$(SRC))
INCLUDES=-I./include #$(addprefix -I,$(SRC_DIR))

vpath %.cpp $(SRC_DIR)

define make-goal
$1/%.o: %.cpp
	$(CC) $(CXXFLAGS) $(INCLUDES) $$< -o $$@
endef

.PHONY: all checkdirs clean

all: checkdirs build

build: $(OBJ)
	$(CC) -shared -Wl,-soname,libcolloid.so -o libcolloid.so $^ -llapack -lblas  -lfftw3

checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR)

$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))
