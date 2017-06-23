#include "../include/skimage.hpp"

skimage::skimage() {
  started = false;
}

skimage::skimage(uint32_t x, uint32_t y, uint32_t f) {
  buffer = new double[x*y*f];
  x_size = x;
  y_size = y;
  f_size = f;
  fxs = x*f;
  totalsize = x*f*y;
  started = true;
  for (auto i=0;i<x_size;i++) {
    for (auto j=0;j<y_size;j++) {
      for (auto k=0;k<f_size;k++) {
        buffer[coords(i,j,k)] = 0;
      }
    }
  }
}

void skimage::init(uint32_t x, uint32_t y, uint32_t f) {
  buffer = new double[x*y*f];
  x_size = x;
  y_size = y;
  f_size = f;
  fxs = x*f;
  totalsize = x*f*y;
  started = true;
  for (auto i=0;i<x_size;i++) {
    for (auto j=0;j<y_size;j++) {
      for (auto k=0;k<f_size;k++) {
        buffer[coords(i,j,k)] = 0;
      }
    }
  }
}

uint32_t skimage::getxsize() { return x_size; }
uint32_t skimage::getysize() { return y_size; }
uint32_t skimage::getfsize() { return f_size; }

double skimage::get(uint32_t x, uint32_t y, uint32_t f) {
  return buffer[coords(x,y,f)];
}

double skimage::get(uint32_t i) {
  return buffer[i];
}

void skimage::get(uint32_t *x, uint32_t *y, uint32_t *f, double *output, uint32_t n) {
  for (auto i=0;i<n;i++) {
    output[i] = buffer[coords(x[i],y[i],f[i])];
  }
}

void skimage::set(uint32_t x, uint32_t y, uint32_t f, double val) {
  buffer[coords(x,y,f)] = val;
}
void skimage::set(uint32_t *x, uint32_t *y, uint32_t *f, double *val, uint32_t n) {
  for (auto i=0;i<n;i++) {
    buffer[coords(x[i],y[i],f[i])] = val[i];
  }
}

double skimage::comp(uint32_t x, uint32_t y, uint32_t f, double val) {
  return val-buffer[coords(x,y,f)];
}
void skimage::comp(uint32_t *x, uint32_t *y, uint32_t *f, double *val, double *output, uint32_t n) {
  for (auto i=0;i<n;i++) {
    output[i] = val[i] - buffer[coords(x[i],y[i],f[i])];
  }
}

double skimage::mean() {
  double runningtotal = 0.0;
  for (auto i=0;i<totalsize;i++) {
    runningtotal += buffer[i];
  }
  return runningtotal/(double)totalsize;
}

double skimage::median() {
  vector <double> data;

  for (auto i=0;i<totalsize;i++) {
    data.push_back(buffer[i]);
  }

  sort(data.begin(),data.end());

  return data[totalsize/2];
}

double skimage::max() {
  double max = 0.0;
  for (auto i=0;i<totalsize;i++) {
    if (buffer[i]>max) { max=buffer[i]; }
  }
  return max;
}

double skimage::min() {
  double min = max();
  for (auto i=0;i<totalsize;i++) {
    if (buffer[i]<min) { min=buffer[i]; }
  }
  return min;
}

void skimage::peak(uint32_t *x, uint32_t *y) {
  double max = 0.0;
  for (auto v=0;v<y_size;v++) {
    for (auto u=0;u<x_size;u++) {
      if (buffer[coords(u,v,0)]>max) {
        max=buffer[coords(u,v,0)];
        *x = u;
        *y = v;
      }
    }
  }
}

void skimage::unnan() {
  for (auto i=0;i<totalsize;i++) {
    if (buffer[i]!=buffer[i]) {
      buffer[i] = 0;
    }
  }
}

double skimage::stdev(double mu) {
  double runningtotal = 0.0;
  for (auto i=0;i<totalsize;i++) {
    runningtotal += (buffer[i]-mu)*(buffer[i]-mu);
  }
  runningtotal /= ((double)totalsize - 1.);
  return sqrt(runningtotal);
}

double skimage::noise(skimage &pbeam) {
  double mu0, cutoff, mu, sig, dev;
  double n;

  mu0 = median();
  cutoff = stdev(mu0)*2;
  cutoff *= cutoff;

  mu = 0;
  n = 0;
  for (auto i=0;i<totalsize;i++) {
    dev = buffer[i]-mu0;
    if ((dev*dev)<cutoff) {
      mu += buffer[i]*pbeam.get(i);
      n+=pbeam.get(i);
    }
  }
  mu /= n;

  sig = 0;
  for (auto i=0;i<totalsize;i++) {
    dev = buffer[i]-mu0;
    if ((dev*dev)<cutoff) {
      sig += pbeam.get(i)*(buffer[i]-mu)*(buffer[i]-mu);
    }
  }
  sig /= (n-1.);

  noiseLevel = sqrt(sig);

  return noiseLevel;
}

double skimage::noise() {
  double mu0, cutoff, mu, sig, dev;
  uint32_t n;

  mu0 = median();
  cutoff = stdev(mu0)*2;
  cutoff *= cutoff;

  mu = 0;
  n = 0;
  for (auto i=0;i<totalsize;i++) {
    dev = buffer[i]-mu0;
    if ((dev*dev)<cutoff) {
      mu += buffer[i];
      n++;
    }
  }
  mu /= (double)n;

  sig = 0;
  for (auto i=0;i<totalsize;i++) {
    dev = buffer[i]-mu0;
    if ((dev*dev)<cutoff) {
      sig += (buffer[i]-mu)*(buffer[i]-mu);
    }
  }
  sig /= ((double)n-1.);

  noiseLevel = sqrt(sig);

  return noiseLevel;
}

void skimage::scan(uint32_t psize, double threshold, double stdev) {
  threshold *= psize*psize;
  stdev *= psize*psize;
  fstream logfile;

  logfile.open("imagescan.txt", fstream::out);
  for (auto y=0;y<y_size;y++) {
    for (auto x=0;x<x_size;x++) {
      logfile << " " << buffer[coords(x,y,0)];
    }
    logfile << endl;
  }

  for (auto y=y_size-psize;y>0;y-=psize) {
    for (auto x=0;x<x_size;x+=psize) {
      double pixel = 0.0;
      for (auto u=0;u<psize;u++) {
        for (auto v=0;v<psize;v++) {
          pixel += buffer[coords(y+v,x+u,0)];
        }
      }
      if (pixel>(threshold+3*stdev)) {
        cout << " *";
      } else {
        if (pixel>(threshold+stdev)) {
          cout << " +";
        } else {
          if (pixel>(threshold+0.25*stdev)) {
            cout << " .";
          } else {
            cout << "  ";
          }
        }
      }
    }
    cout << endl;
  }
}

//Crops the image to the centre part
void skimage::crop(uint32_t nx, uint32_t ny) {
  uint32_t dx = (x_size-nx)/2;
  uint32_t dy = (y_size-ny)/2;
  double *newbuffer = new double[nx*ny*f_size];
  double *oldbuffer = buffer;

  for (auto f=0;f<f_size;f++) {
    for (auto u=0;u<nx;u++) {
      for (auto v=0;v<ny;v++) {
        uint32_t position = f*u*f_size+v*nx*f_size;
        uint32_t source = f*(u+dx)*f_size+(v+dy)*ny*f_size;
        newbuffer[position] = oldbuffer[source];
      }
    }
  }

  buffer = newbuffer;
  x_size = nx;
  y_size = ny;

  totalsize = x_size*y_size*f_size;
  fxs = x_size*f_size;

  delete oldbuffer;
}

//Use to pad out an image to avoid having to have boundary checks in map
void skimage::pad(uint32_t dx, uint32_t dy, double contents) {
  uint32_t nx = x_size+dx;
  uint32_t ny = y_size+dy;
  double *newbuffer = new double[nx*ny*f_size];
  double *oldbuffer = buffer;

  dx*=0.5;
  dy*=0.5;

  for (auto f=0;f<f_size;f++) {
    for (auto u=0;u<nx;u++) {
      for (auto v=0;v<ny;v++) {
        uint32_t position = f+u*f_size+v*nx*f_size;
        if (u<dx) { newbuffer[position]=contents; continue; }
        if (u>(nx-dx-1)) { newbuffer[position]=contents; continue; }
        if (v<dy) { newbuffer[position]=contents; continue; }
        if (v>(ny-dy-1)) { newbuffer[position]=contents; continue; }
        newbuffer[position] = oldbuffer[coords(u-dx,v-dy,f)];
      }
    }
  }

  buffer = newbuffer;
  x_size = nx;
  y_size = ny;

  totalsize = x_size*y_size*f_size;
  fxs = x_size*f_size;

  delete oldbuffer;
}

double skimage::deconv(skimage &other, twov *points, double *flux, uint32_t n) {
  double result = 0;

  uint32_t cx = other.getxsize()*0.5;
  uint32_t cy = other.getysize()*0.5;
  for (auto i=0;i<n;i++) {
    result += 2*flux[i]*buffer[coords(points[i].x, points[i].y,0)];
    for (auto j=0;j<n;j++) {
      uint32_t u = cx + points[j].x-points[i].x;
      uint32_t v = cy + points[j].y-points[i].y;
      result -= flux[i]*flux[j]*other.get(u,v,0);
    }
  }

  return result * 0.5/(noiseLevel*noiseLevel);
}

double skimage::deconv(skimage &other, skimage &pbeam, twov *points, double *flux, uint32_t n) {
  double result = 0;
  double ourflux[n];

  for (auto i=0;i<n;i++) {
    ourflux[i] = flux[i]*pbeam.get(points[i].x,points[i].y,0);
  }

  uint32_t cx = other.getxsize()*0.5;
  uint32_t cy = other.getysize()*0.5;
  for (auto i=0;i<n;i++) {
    result += 2*ourflux[i]*buffer[coords(points[i].x, points[i].y,0)];
    for (auto j=0;j<n;j++) {
      uint32_t u = cx + points[j].x-points[i].x;
      uint32_t v = cy + points[j].y-points[i].y;
      result -= ourflux[i]*ourflux[j]*other.get(u,v,0);
    }
  }

  return result * 0.5/(noiseLevel*noiseLevel);
}

void skimage::subbeam(skimage &other, twov points, double flux) {
  uint32_t cx = other.getxsize()*0.5;
  uint32_t cy = other.getysize()*0.5;

  for (auto y=0;y<y_size;y++) {
    for (auto x=0;x<x_size;x++) {
        uint32_t u = cx + points.x-x;
        uint32_t v = cy + points.y-y;
        buffer[coords(x,y,0)] -= other.get(u,v,0);
    }
  }
}

void skimage::subtract(twov points, double flux) {
  buffer[coords(points.x,points.y,0)] -= flux;
}

void skimage::add(twov points, double flux) {
  buffer[coords(points.x,points.y,0)] += flux;
}

// This was to get an estimate, this method is not correct
double skimage::badresidual(double mean, double noise) {
  double result = 0;

  for (auto y=0;y<y_size;y++) {
    for (auto x=0;x<x_size;x++) {
      double resid = (buffer[coords(x,y,0)] - mean)/noise;
      result += resid*resid;
    }
  }

  return sqrt(result/(x_size*y_size));
}

double skimage::fluxscale(skimage &beam) {
  //cout << max() << " " << beam.max() << endl;
  return 2*max()/beam.max();
}

void skimage::setnoise(double value) {
  noiseLevel = value;
}

double skimage::getnoise() {
  return noiseLevel;
}
