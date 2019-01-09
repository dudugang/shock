CXX = g++
exe = flow
src = src
obj = obj
bin = bin
build = build
include_paths = include include/eigen
include = $(foreach dir, $(include_paths), -I$(dir))

.PHONY: all
all: directories $(bin)/$(exe)

# This is purely for testing purposes
.PHONY: print
print:
	$(info $(include))

.PHONY: directories
directories:
	if [ ! -d bin ]; then mkdir bin; fi
	if [ ! -d obj ]; then mkdir obj; fi

.PHONY: clean
clean:
	rm -f *.o $(obj)/* $(bin)/$(exe)
	rmdir obj/ bin/

$(bin)/$(exe): $(obj)/Algebra.o $(obj)/Flux.o $(obj)/Flow.o $(obj)/main.o
	$(CXX) -o $@ $(include) $(obj)/Algebra.o $(obj)/Flux.o $(obj)/Flow.o $(obj)/main.o

$(obj)/Algebra.o: $(src)/Algebra.cpp
	$(CXX) -o $@ $(include) -c $<

$(obj)/Flux.o: $(src)/Flux.cpp
	$(CXX) -o $@ $(include) -c $<

$(obj)/Flow.o: $(src)/Flow.cpp
	$(CXX) -o $@ $(include) -c $<

$(obj)/main.o: $(src)/main.cpp
	$(CXX) -o $@ $(include) -c $<
