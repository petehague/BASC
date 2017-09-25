/*
  Options handler.

  Works with both option files and command line arguments
  Assumes format KEY=value
*/

#include <vector>
#include <string>
#include <cstring>
#include <cinttypes>
#include <fstream>
#include <iostream>

#include "../include/options.hpp"

using namespace std;

/*
  Parses a single item KEY=value
  if it encounters a single work with no '=' this is assumed to be
  an option filename, and the file is opened and parsed
*/
void optDict::parse(string item) {
    char c_item[item.size()];
    string *valstr;
  	char *keyname, *value;

    strcpy(c_item, item.c_str());
	  keyname = strtok(c_item, "= ");
	  if (keyname!=NULL && keyname[0]!='#') {
        name.push_back(keyname);
        value = strtok(NULL, "= ");
        if (value == NULL) {
            readFile(keyname);
        } else {
        		valstr = new string(value);
            if ((*valstr).find_first_not_of("0123456789+-.")==string::npos) {
                if ((*valstr).find_first_of(".")==string::npos) {
                    data.push_back(new int32_t(stoi(*valstr)));
                    type.push_back(opt_int);
                    delete valstr;
                } else {
                    data.push_back(new double(stod(*valstr)));
                    type.push_back(opt_real);
                    delete valstr;
                }
            } else {
                data.push_back((void *)valstr);
                type.push_back(opt_str);
            }
        }
    }
}

void optDict::readFile (string filename) {
    fstream optfile;

    optfile.open(filename.c_str(), fstream::in);
    for (string line; getline(optfile, line); ) {
        parse(line);
    }
    optfile.close();
}

void optDict::readCL (int argc, char **argv){
  	for (int i=1;i<argc;i++) {
        parse(argv[i]);
    }
}

/*
  Retrival routines - these rely on the code getting the type correct.
  May add a check in later if this becomes a problem
*/
int32_t optDict::getint(string key) {
    int index = 0;
    while (key.compare(name[index])!=0 && index<data.size()) index++;
    if (index==data.size()) return 0;
    return *(int32_t *)data[index];
}

double optDict::getdouble(string key) {
    int index = 0;
    while (key.compare(name[index])!=0 && index<data.size()) index++;
    if (index==data.size()) return 0;
    return *(double *)data[index];
}

string optDict::getstring(string key) {
    int index = 0;
    while (key.compare(name[index])!=0 && index<data.size()) index++;
    if (index==data.size()) return "";
    return *(string *)data[index];
}

void optDict::report() {
    for (uint32_t i=0;i<name.size();i++) {
        cout << name[i] << " = ";
        switch (type[i]) {
            case opt_str:
                cout << *(string *)data[i];
                break;

            case opt_int:
                cout << *(int32_t *)data[i];
                break;

            case opt_real:
                cout << *(double *)data[i];
                break;
        }
        cout << endl;
    }
}
