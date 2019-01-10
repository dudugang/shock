CXX = g++
exe = flow
src_dir = src
src = $(wildcard $(src_dir)/*.cpp)
bin_dir = bin
build_dir = build
obj_dir = $(build_dir)/obj
obj = $(src:$(src_dir)/%.cpp=$(obj_dir)/%.o)
test_dir = test
test_obj_dir = $(build_dir)/test
include_paths = include lib/eigen lib/googletest
cpp_flags = $(foreach dir, $(include_paths), -I$(dir))
c_flags = -O3 -Wall

.PHONY: all
all: directories $(bin_dir)/$(exe)

# This is purely for testing purposes
.PHONY: print
print:
	$(info $(obj))

.PHONY: directories
directories:
	if [ ! -d $(bin_dir) ]; then mkdir $(bin_dir); fi
	if [ ! -d $(build_dir) ]; then mkdir $(build_dir); fi
	if [ ! -d $(obj_dir) ]; then mkdir $(obj_dir); fi
	if [ ! -d $(test_obj_dir) ]; then mkdir $(test_obj_dir); fi

.PHONY: clean
clean:
	rm -rf $(build_dir) $(bin_dir)

$(bin_dir)/$(exe): $(obj)
	$(CXX) -o $@ $^

$(obj_dir)/%.o: $(src_dir)/%.cpp
	$(CXX) $(cpp_flags) $(c_flags) -o $@ -c $<
