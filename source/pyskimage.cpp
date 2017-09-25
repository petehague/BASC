#include <Python.h>

#include "../include/skimage.hpp"

skimage dirtyMap, dirtyBeam;

//Python stuff

twov getvec(PyObject *tup) {
  twov result;
  double flux;
  uint32_t x,y;
  PyArg_ParseTuple(tup,"iid", &x, &y, &flux);
  result.x = (double)x;
  result.y = (double)y;
  return result;
}

double getflux(PyObject *tup) {
  double flux;
  uint32_t x,y;
  PyArg_ParseTuple(tup,"iid", &x, &y, &flux);
  return flux;
}


void parsefrompy(skimage *target, PyObject* source, uint32_t x, uint32_t y) {
  PyObject *datum;
  uint32_t counter;
  double value;

  target->init(x,y,1);

  if (PyList_Size(source)<(x*y)) {
    cout << "Not enough data provided" << endl;
    return;
  }

  for (uint32_t v=0;v<y;v++) {
    for (uint32_t u=0;u<x;u++) {
      counter = u+v*x;
      datum = PyList_GetItem(source, counter);
      PyArg_Parse(datum, "d", &value);
      target->set(u,v,0,value);
    }
  }
}

static PyObject* skimage_makebeam(PyObject *self, PyObject *args) {
  uint32_t xb,yb;
  PyObject *beam;

  if (!PyArg_ParseTuple(args,"O!ii", &PyList_Type, &beam, &xb, &yb)) { return NULL; }

  parsefrompy(&dirtyBeam,beam,xb,yb);

  dirtyBeam.pad(xb,yb,0);

  return PyInt_FromLong(1);
}

static PyObject* skimage_makemap(PyObject *self, PyObject *args) {
  uint32_t xm,ym;
  PyObject *map;

  if (!PyArg_ParseTuple(args,"O!ii", &PyList_Type, &map, &xm, &ym)) { return NULL; }

  parsefrompy(&dirtyMap,map,xm,ym);

  return PyInt_FromLong(1);
}

static PyObject* skimage_deconvolve(PyObject *self, PyObject *args) {
  PyObject *points;
  PyObject *datum;
  uint32_t psize;
  vector <twov> atoms;
  vector <double> flux;
  double result = 0;

  if (!PyArg_ParseTuple(args,"O!", &PyList_Type, &points)) { return NULL; }

  psize = (uint32_t)PyList_Size(points);
  for (uint32_t i=0;i<psize;i++) {
    datum = PyList_GetItem(points,i);
    atoms.push_back(getvec(datum));
    flux.push_back(getflux(datum));
  }

  result = dirtyMap.deconv(dirtyBeam, &atoms[0], &flux[0], psize);

  return Py_BuildValue("d",result);
}

static PyObject* skimage_noise(PyObject *self, PyObject *args) {
  //double mu = dirtyMap.mean();
  double stdev;

  //stdev = dirtyMap.stdev(mu);
  stdev = dirtyMap.noise();

  return Py_BuildValue("d", stdev);
}

static PyObject* skimage_render(PyObject *self, PyObject *args) {
  uint32_t x,y;
  skimage image;
  PyObject *data;

  if (!PyArg_ParseTuple(args,"O!ii", &PyList_Type, &data, &x, &y)) { return NULL; }

  parsefrompy(&image, data, x, y);

  image.scan(16, image.mean(),image.stdev(image.mean()));

  return PyInt_FromLong(1);
}

static PyObject* skimage_fluxscale(PyObject *self, PyObject *args) {
  dirtyBeam.scan(32, dirtyBeam.mean(), dirtyBeam.stdev(dirtyBeam.mean()));
  return Py_BuildValue("d", dirtyMap.fluxscale(dirtyBeam));
}

static PyObject* skimage_version(PyObject *self, PyObject *args) {
  cout << "Version number 0.1" << endl;
  return PyLong_FromLong(1);
}

static PyMethodDef skmodule[] = {
  {"version", skimage_version, METH_VARARGS, ""},
  {"render", skimage_render, METH_VARARGS, ""},
  {"deconvolve", skimage_deconvolve, METH_VARARGS, ""},
  {"makeMap", skimage_makemap, METH_VARARGS, ""},
  {"makeBeam", skimage_makebeam, METH_VARARGS, ""},
  {"mapNoise", skimage_noise, METH_VARARGS,""},
  {"fluxScale", skimage_fluxscale, METH_VARARGS,""},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initskimage(void) {
#ifdef PYTHREE
  return PyModule_Create(&skmodule);
#else
  (void) Py_InitModule("skimage", skmodule);
#endif
}
