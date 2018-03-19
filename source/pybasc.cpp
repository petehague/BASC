#ifndef STANDALONE
#include <Python.h>
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <vector>

#include "../include/skimage.hpp"
#include "../bayesys/bayesys3.h"
#include "../include/options.hpp"

using namespace std;

vector <skimage> dirtyMap;
vector <skimage> dirtyBeam;
vector <skimage> primaryBeam;
vector <double> evidence;
uint32_t mapsize, mapdepth, modelIndex, maxmodels;
double sigma, freqscale;
bool cubeSwitch = false;
optDict options;

struct UserCommonStr {
  uint32_t nmodels;
  uint32_t maxmodels;
  uint32_t burnin;
  uint32_t mapindex;
  double noise;
  double fluxscale;
  double cool;
  vector <uint16_t> natoms;
  vector <double> x;
  vector <double> y;
  vector <double> F;
  vector <double> fmu;
  vector <double> fsig;
  vector <double> lh;
  vector <uint32_t> modelIndex;
};

vector <UserCommonStr> UserCommon;
CommonStr Common;

vector <double> fx;
vector <double> fF;
vector <double> fL;

extern "C" {
int UserBuild(double *like, CommonStr *Common, ObjectStr* Member, int natoms, int dummy) {
  vector <twov> points;
  vector <double> flux;
  vector <uint32_t> ax, ay;
  vector <double> fmu;
  vector <double> fsig;
  double **Cube = Member->Cubes;
  UserCommonStr *UC = (UserCommonStr *)Common->UserCommon;
  double fluxscale = UC->fluxscale;

  if (natoms==0) {
    *like = -1e6;
    return 0;
    //*like = 100;
    //return 1;
  }

  for (uint32_t i=0;i<natoms;i++) {
    double x = Cube[i][0]*double(mapsize);
    double y = Cube[i][1]*double(mapsize);
    points.push_back(twov(x,y));
    ax.push_back(Cube[i][0]*mapsize);
    ay.push_back(Cube[i][1]*mapsize);
    flux.push_back(fluxscale*Cube[i][2]/(1.-Cube[i][2]));
    //fx.push_back(Cube[i][2]);
    //fF.push_back(flux.back());
    if (cubeSwitch) {
      fmu.push_back(Cube[i][3]*freqscale);
      fsig.push_back(Cube[i][4]*freqscale);
    }
  }

  if (cubeSwitch) {
    *like = dirtyMap[UC->mapindex].deconv(dirtyBeam[UC->mapindex],primaryBeam[UC->mapindex],&points[0],&flux[0],&fmu[0],&fsig[0],flux.size());
    return 1;
  }

  *like = dirtyMap[UC->mapindex].deconv(dirtyBeam[UC->mapindex],primaryBeam[UC->mapindex],&ax[0], &ay[0],&flux[0],flux.size());
  //*like = dirtyMap[UC->mapindex].deconv(dirtyBeam[UC->mapindex],&points[0],&flux[0],flux.size());

  /*for (uint32_t i=0;i<natoms;i++) {
    fL.push_back(*like);
  }*/

  return 1;
}

int UserMonitor(CommonStr *Common, ObjectStr *Members) {
  double **Cube;// = Members->Cubes;
  UserCommonStr *UC = (UserCommonStr *)Common->UserCommon;
  double fluxscale = UC->fluxscale;

  if (Common->cool > 1.) Common->cool = 1.;
  if (Common->cool < 1.) {
    UC->burnin += 1;
    if (Common->cool>UC->cool) {
      //cout << Common->cool << " " << UC->burnin << endl;
      UC->cool += 0.1;
    }
    return 0;
  }

  if (UC->maxmodels>UC->burnin) { UC->maxmodels = UC->burnin+1; }

  UC->natoms.push_back(Members[0].Natoms);
  if (UC->nmodels == 0) {
    cout << "Burn in complete" << endl;
  }
  UC->nmodels += 1;
  for (uint32_t k=0;k<Common->ENSEMBLE;k++) {
    Cube = Members[k].Cubes;
    for (uint32_t i=0;i<Members[k].Natoms;i++) {
      UC->lh.push_back(Members[k].Lhood);
      UC->x.push_back(Cube[i][0]*mapsize);
      UC->y.push_back(Cube[i][1]*mapsize);
      UC->F.push_back(fluxscale*(Cube[i][2]/(1.-Cube[i][2])));
      if (cubeSwitch) {
        UC->fmu.push_back(Cube[i][3]*freqscale);
        UC->fsig.push_back(Cube[i][4]*freqscale);
      }
      UC->modelIndex.push_back(modelIndex);
    }
    modelIndex++;
  }

  if (UC->nmodels>UC->maxmodels) {
    cout << "Chain complete with " << UC->nmodels << " models" << endl;
    return 1;
  } else {
    return 0;
  }
}
}

#ifndef STANDALONE
void parsefrompy(skimage &target, PyObject* source, uint32_t x, uint32_t y) {
  PyObject *datum;
  uint32_t counter;
  double value;

  target.init(x,y,1);

  if (PyList_Size(source)<(x*y)) {
    cout << "Not enough data provided" << endl;
    return;
  }

  for (uint32_t v=0;v<y;v++) {
    for (uint32_t u=0;u<x;u++) {
      counter = u+v*x;
      datum = PyList_GetItem(source, counter);
      PyArg_Parse(datum, "d", &value);
      target.set(u,v,0,value);
    }
  }
}


static PyObject* bascmodule_setgrid(PyObject *self, PyObject *args) {
  uint32_t crpix1, crpix2, mapid, cindex;
  double crval1, crval2, cdelt1, cdelt2;

  if (!PyArg_ParseTuple(args,"iiiddidd", &cindex, &mapid, &crpix1, &crval1, &cdelt1, &crpix2, &crval2, &cdelt2)) { return NULL; }

  if (mapid==0) {
    dirtyMap[cindex].setgrid(0,crpix1,crval1,cdelt1);
    dirtyMap[cindex].setgrid(1,crpix2,crval2,cdelt2);
  }

  if (mapid==1) {
    dirtyBeam[cindex].setgrid(0,crpix1,crval1,cdelt1);
    dirtyBeam[cindex].setgrid(1,crpix2,crval2,cdelt2);
  }

  if (mapid==2) {
    primaryBeam[cindex].setgrid(0,crpix1,crval1,cdelt1);
    primaryBeam[cindex].setgrid(1,crpix2,crval2,cdelt2);
  }

  return PyLong_FromLong(1);
}

static PyObject* bascmodule_map(PyObject *self, PyObject *args) {
  uint32_t xb,yb,id,cindex;
  PyObject *map;

  if (!PyArg_ParseTuple(args,"iO!iii", &cindex, &PyList_Type, &map, &xb, &yb, &id)) { return NULL; }

  if (id==0) {
    parsefrompy(dirtyMap.at(cindex),map,xb,yb);
    dirtyMap[cindex].crop(xb/2, yb/2);
  }

  if (id==1) {
    parsefrompy(dirtyBeam.at(cindex),map,xb,yb);
  }

  if (id==2) {
    parsefrompy(primaryBeam.at(cindex),map,xb,yb);
    primaryBeam[cindex].crop(xb/2, yb/2);
    primaryBeam[cindex].unnan();
  }

  //TODO: add support for non square maps
  mapsize = xb/2;

  //TODO: add support for Cubes
  mapdepth = 1;

  return PyLong_FromLong(1);
}

static PyObject* bascmodule_getmap(PyObject *self, PyObject *args) {
  uint32_t mapIndex, listSize, axisSize, cindex;
  PyObject *returnList;
  skimage mapRef;

  if (!PyArg_ParseTuple(args,"ii", &cindex, &mapIndex)) { return NULL; }

  if (mapIndex==1) {
    axisSize = mapsize*2;
    mapRef = dirtyBeam.at(cindex);
  } else {
    axisSize = mapsize;
    if (mapIndex==0) { mapRef = dirtyMap.at(cindex); }
    if (mapIndex==2) { mapRef = primaryBeam.at(cindex); }
  }

  listSize = axisSize*axisSize;
  returnList = PyList_New(listSize);
  for (uint32_t y=0;y<axisSize;y++) {
      for (uint32_t x=0;x<axisSize;x++) {
        PyList_SetItem(returnList, x+y*axisSize, PyFloat_FromDouble(mapRef.get(x,y,0)));
      }
  }

  return returnList;
}

static PyObject* bascmodule_init(PyObject *self, PyObject *args) {
  //Common.UserCommon = (void *)UserCommon;

  return PyLong_FromLong(1);
}

static PyObject* bascmodule_setOption(PyObject *self, PyObject *args) {
  char *key;
  char *value;

  if (!PyArg_ParseTuple(args,"ss", &key, &value)) { return NULL; }

  if (strcmp(key,"MinAtoms")==0) { Common.MinAtoms = atoi(value); }
  if (strcmp(key,"MaxAtoms")==0) { Common.MaxAtoms = atoi(value); }
  if (strcmp(key,"Alpha")==0) { Common.Alpha = atoi(value); }
  if (strcmp(key,"Valency")==0) { Common.Valency = atoi(value); }
  if (strcmp(key,"Iseed")==0) { Common.Iseed = atoi(value); }
  if (strcmp(key,"Ensemble")==0) { Common.ENSEMBLE = atoi(value); }
  if (strcmp(key,"Method")==0) { Common.Method = atoi(value); }
  if (strcmp(key,"Rate")==0) { Common.Rate = atof(value); }
  if (strcmp(key,"maxmodels")==0) { maxmodels = atoi(value); }

  return PyLong_FromLong(1);
}

static PyObject* bascmodule_setNoise(PyObject *self, PyObject *args) {
  uint32_t index;
  double newnoise;

  if (!PyArg_ParseTuple(args,"id", &index, &newnoise)) { return NULL; }

  UserCommon[index].noise = newnoise;

  return PyLong_FromLong(1);
}

static PyObject* bascmodule_setFlux(PyObject *self, PyObject *args) {
  uint32_t index;
  double newflux;

  if (!PyArg_ParseTuple(args,"id", &index, &newflux)) { return NULL; }

  UserCommon[index].fluxscale = newflux;

  return PyLong_FromLong(1);
}

static PyObject* bascmodule_run(PyObject *self, PyObject *args) {
  int32_t code, cindex;
  ObjectStr *Members;
  fstream chainfile;
  double noiseflux;

  if (!PyArg_ParseTuple(args,"i", &cindex)) { return NULL; }

  noiseflux = dirtyMap[cindex].noise(primaryBeam.at(cindex));
  dirtyMap[cindex].setnoise(noiseflux);

  if (UserCommon[cindex].noise>0) {
     dirtyMap[cindex].setnoise(UserCommon[cindex].noise);
     noiseflux = UserCommon[cindex].noise;
  }
  if (UserCommon[cindex].fluxscale==0) {
     UserCommon[cindex].fluxscale = noiseflux;
     cout << "Defaulting to flux scaled to noise." << endl;
  } 

  //cout << "PSF min/max " << dirtyBeam[cindex].min() << "," << dirtyBeam[cindex].max() << endl;

  cout << "Using noise level " << dirtyMap[cindex].getnoise() << endl;
  cout << "Flux prior centre " << UserCommon[cindex].fluxscale << endl;

  //dirtyMap[cindex].scan(64,0,0.001);
  //dirtyBeam[cindex].scan(128, 0, 0.01);

  //cout << "Dirty map maximum " << dirtyMap[cindex].max() << endl;

  if (mapdepth>1) {
    Common.Ndim = 5;
  } else {
    Common.Ndim = 3;
  }

  Common.cool = 0;

  UserCommon[cindex].nmodels = 0;
  UserCommon[cindex].burnin = 0;
  UserCommon[cindex].cool = 0;
  UserCommon[cindex].maxmodels = maxmodels;
  UserCommon[cindex].mapindex = cindex;
  Common.UserCommon = (void *)&(UserCommon.data()[cindex]);

  Members = new ObjectStr[Common.ENSEMBLE];
  code = BayeSys3(&Common,Members);

  evidence[cindex] = Common.Evidence;

  /*chainfile.open("chain.txt", fstream::out);
  for (uint32_t i=0;i<UserCommon[cindex].x.size();i++) {
    chainfile << fixed << setprecision(9) << UserCommon[cindex].x[i] << " " << UserCommon[cindex].y[i] << " " << UserCommon[cindex].F[i] << " " << UserCommon[cindex].modelIndex[i];
    chainfile << endl;
  }
  chainfile.close();*/

  return PyLong_FromLong(code);
}

static PyObject* bascmodule_evid(PyObject *self, PyObject *args) {
  uint32_t cindex;

  if (!PyArg_ParseTuple(args,"i", &cindex)) { return NULL; }

  return PyFloat_FromDouble(evidence[cindex]);
}

static PyObject* bascmodule_chain(PyObject *self, PyObject *args) {
  uint32_t colIndex, cindex;
  PyObject *returnList;
  vector <double> colPtr;

  if (!PyArg_ParseTuple(args,"ii", &cindex, &colIndex)) { return NULL; }
  returnList = PyList_New(UserCommon[cindex].nmodels);

  if (colIndex==0) { colPtr=UserCommon[cindex].x; }
  if (colIndex==1) { colPtr=UserCommon[cindex].y; }
  if (colIndex==2) { colPtr=UserCommon[cindex].F; }
  if (colIndex==3) {
    vector <uint32_t> intPtr=UserCommon[cindex].modelIndex;
    for (uint32_t i=0;i<UserCommon[cindex].nmodels;i++) {
      PyList_SetItem(returnList, i, PyLong_FromLong(intPtr[i]));
    }
    return returnList;
  }
  if (colIndex==4) { colPtr=UserCommon[cindex].lh; }

  for (uint32_t i=0;i<UserCommon[cindex].nmodels;i++) {
    PyList_SetItem(returnList, i, PyFloat_FromDouble(colPtr[i]));
  }

  return returnList;
}

static PyObject* bascmodule_show(PyObject *self, PyObject *args) {
  uint32_t item, cindex;
  skimage skRef;

  if (!PyArg_ParseTuple(args,"ii", &cindex, &item)) { return NULL; }

  if (item==0) { skRef = dirtyMap.at(cindex); }
  if (item==1) { skRef = dirtyBeam.at(cindex); }
  if (item==2) { skRef = primaryBeam.at(cindex); }

  skRef.scan(32, skRef.mean(), 0.5*skRef.stdev(skRef.mean()));
  return PyLong_FromLong(1);
}

static PyObject* bascmodule_version(PyObject *self, PyObject *args) {
  cout << "Version number 0.1" << endl;
  return PyLong_FromLong(1);
}

static PyObject* bascmodule_new(PyObject *self, PyObject *args) {
  dirtyMap.push_back(skimage());
  dirtyBeam.push_back(skimage());
  primaryBeam.push_back(skimage());
  UserCommon.push_back(UserCommonStr());
  UserCommon.back().noise = -1;
  UserCommon.back().fluxscale = 0;
  evidence.push_back(-1.0);
  return PyLong_FromLong(dirtyMap.size()-1);
}

static PyMethodDef bascmethods[] = {
  {"version", bascmodule_version, METH_VARARGS, ""},
  {"grid", bascmodule_setgrid, METH_VARARGS, ""},
  {"map", bascmodule_map, METH_VARARGS, ""},
  {"getmap", bascmodule_getmap, METH_VARARGS, ""},
  {"run", bascmodule_run, METH_VARARGS, ""},
  {"chain", bascmodule_chain, METH_VARARGS, ""},
  {"show", bascmodule_show, METH_VARARGS, ""},
  {"init", bascmodule_init, METH_VARARGS, ""},
  {"option", bascmodule_setOption, METH_VARARGS, ""},
  {"noise", bascmodule_setNoise, METH_VARARGS, ""},
  {"flux", bascmodule_setFlux, METH_VARARGS, ""},
  {"new", bascmodule_new, METH_VARARGS, ""},
  {"evidence", bascmodule_evid, METH_VARARGS, ""},
  {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION > 2
static struct PyModuleDef bascmodule = {
  PyModuleDef_HEAD_INIT,
  "bascmod",
  NULL,
  -1,
  bascmethods
};
PyMODINIT_FUNC
PyInit_bascmod(void) {
  return PyModule_Create(&bascmodule);
#else
PyMODINIT_FUNC
initbascmod(void) {
  (void) Py_InitModule("bascmod", bascmethods);
#endif
}
#else
int main() {
  CommonStr Common;
  ObjectStr *Members;

  dirtyMap.push_back(skimage());
  dirtyBeam.push_back(skimage());
  primaryBeam.push_back(skimage());
  UserCommon.push_back(UserCommonStr());

  dirtyMap[0].init(10,10,1);
  dirtyBeam[0].init(20,20,1);
  primaryBeam[0].init(10,10,1);

  options.readFile("config.txt");

  Common.Ndim = 3;
  Common.MinAtoms = options.getint("MinAtoms");
  Common.MaxAtoms = options.getint("MaxAtoms");
  Common.Alpha = options.getint("Alpha");
  Common.Valency = options.getint("Valency");
  Common.Iseed =  options.getint("Iseed");
  Common.ENSEMBLE =  options.getint("Ensemble");
  Common.Method =  options.getint("Method");
  Common.Rate =  options.getint("Rate");
  Common.UserCommon = (void *)&(UserCommon.at(0));
  UserCommon[0].nmodels = 0;
  UserCommon[0].maxmodels =  options.getint("maxmodels");
  UserCommon[0].burnin = 0;
  UserCommon[0].cool = 0;
  UserCommon[0].fluxscale = 0.0025;

  Members = new ObjectStr[Common.ENSEMBLE];

  cout << BayeSys3(&Common,Members) << endl;
}
#endif
