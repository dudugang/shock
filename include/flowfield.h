#include <vector>
#include <unordered_map>
#include <input.h>
#include <cell.h>
using std::vector;
using std::unordered_map;

class Flowfield {
    public:
        Flowfield(Input);
        Input input;
        unordered_map<int, Cell> cells;
};
