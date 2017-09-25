#ifndef D_SKIMAGE
#define D_SKIMAGE
#include <iostream>
#include <cinttypes>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

struct twov {
  double x;
  double y;
  double dot(twov other) { return x*other.x + y*other.y; }
  void add(twov other) { x+=other.x; y+=other.y; }
  twov() { x=0; y=0; }
  twov(double xx, double yy) { x=xx; y=yy; }
};

class skimage {
  uint32_t x_size, y_size, f_size, fxs, totalsize;
  double x_0, y_0, x_delt, y_delt;
  //double freq_ref, dfreq;
  double *buffer;
  uint32_t coords(uint32_t x, uint32_t y, uint32_t f) {
    return f+(x*f_size)+(y*fxs);
  }
  bool started;
  double noiseLevel;

public:
  skimage();
  skimage(uint32_t x, uint32_t y, uint32_t f);
  void init(uint32_t x, uint32_t y, uint32_t f);
  void unnan();
  uint32_t getxsize();
  uint32_t getysize();
  uint32_t getfsize();
  double get(uint32_t x, uint32_t y, uint32_t f);
  double get(uint32_t i);
  void get(uint32_t *x, uint32_t *y, uint32_t *f, double *output, uint32_t n);
  void set(uint32_t x, uint32_t y, uint32_t f, double val);
  void set(uint32_t *x, uint32_t *y, uint32_t *f, double *val, uint32_t n);
  double comp(uint32_t x, uint32_t y, uint32_t f, double val);
  void comp(uint32_t *x, uint32_t *y, uint32_t *f, double *val, double *output, uint32_t n);
  double mean();
  double median();
  double max();
  double min();
  void peak(uint32_t *x, uint32_t *y);
  double stdev(double mu);
  double noise();
  double noise(skimage &pbeam);
  void scan(uint32_t psize, double threshold, double stdev);
  void crop(uint32_t nx, uint32_t ny);
  void pad(uint32_t dx, uint32_t dy, double contents);
  double deconv(skimage &other, twov *points, double *flux, uint64_t n);
  double deconv(skimage &other, skimage &pbeam, twov *points, double *flux, double *fmu, double *fsig, uint64_t n);
  double deconv(skimage &other, skimage &pbeam, twov *points, double *flux, uint64_t n);
  void subbeam(skimage &other, twov points, double flux);
  void subtract(twov points, double flux);
  void add(twov points, double flux);
  double badresidual(double mean, double noise);
  double fluxscale(skimage &beam);
  void setnoise(double value);
  double getnoise();
  void setgrid(uint8_t axis, uint32_t crpix, double crval, double cdelt);
  double real(uint8_t axis, double value);
};
#endif
