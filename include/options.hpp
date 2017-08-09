#include <vector>
#include <string>
#include <cinttypes>

using namespace std;

enum typeID {
    opt_str,
    opt_real,
    opt_int
};

class optDict {
    vector <string> name;
    vector <void *> data;
    vector <typeID> type;
public:
    optDict () { }
    optDict (string filename) { readFile(filename); }
    void parse(string item);
    void readFile (string filename);
    void readCL (int argc, char **argv);
    int32_t getint(string key);
    double getdouble(string key);
    string getstring(string key);
    void report();
};
