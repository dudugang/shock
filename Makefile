# Compiler
CXX = g++
# Executable names
exe = shock
test_exe = shocktest
# Directories
src_dir = src
bin_dir = bin
build_dir = build
test_dir = test
obj_dir = $(build_dir)/obj
test_obj_dir = $(test_dir)/obj
test_src_dir = $(test_dir)/src
# Files
test_src = $(wildcard $(test_src_dir)/*.cpp)
src = $(wildcard $(src_dir)/*.cpp)
obj = $(src:$(src_dir)/%.cpp=$(obj_dir)/%.o)
test_obj = $(test_src:$(test_src_dir)/%.cpp=$(test_obj_dir)/%.o)
# Paths to includes
include_paths = include lib/googletest/googletest/include lib/hdf5/hdf5/include lib/hdf5/hdf5/lib
# Compiler flags
flags = $(foreach dir, $(include_paths), -I$(dir)) -std=c++11 -g -Wall
# Libraries and locations
ldlibs = -Llib/googletest/lib -Llib/hdf5/hdf5/lib -lgtest -lgtest_main -lpthread -lhdf5 -lhdf5_cpp -Wl,-rpath=lib/hdf5/hdf5/lib
# Useful variables
empty =
test_suffix = _test

.PHONY: all
all: directories $(bin_dir)/$(exe)

.PHONY: test
test: test_directories $(bin_dir)/$(test_exe)

# This is purely for testing purposes
.PHONY: print
print:
	$(info $(patsubst $(test_obj_dir)/%_test.o,$(obj_dir)/%.o,$(test_obj)))

.PHONY: directories
directories:
	if [ ! -d $(bin_dir) ]; then mkdir $(bin_dir); fi
	if [ ! -d $(build_dir) ]; then mkdir $(build_dir); fi
	if [ ! -d $(obj_dir) ]; then mkdir $(obj_dir); fi

.PHONY: test_directories
test_directories:
	if [ ! -d $(test_dir) ]; then mkdir $(test_dir); fi
	if [ ! -d $(test_obj_dir) ]; then mkdir $(test_obj_dir); fi

.PHONY: clean
clean:
	rm -rf $(build_dir) $(bin_dir) $(test_obj_dir)

$(bin_dir)/$(exe): $(obj)
	$(CXX) -o $@ $^ $(ldlibs)

$(obj_dir)/%.o: $(src_dir)/%.cpp
	$(CXX) $(flags) -o $@ -c $<

$(bin_dir)/$(test_exe): $(patsubst $(test_obj_dir)/%_test.o,$(obj_dir)/%.o,$(test_obj)) $(test_obj)
	$(CXX) -o $@ $^ $(ldlibs)

$(test_obj_dir)/%.o: $(test_src_dir)/%.cpp
	$(CXX) $(flags) -o $@ -c $<
