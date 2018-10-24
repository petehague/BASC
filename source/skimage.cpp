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
  for (uint32_t i=0;i<x_size;i++) {
    for (uint32_t j=0;j<y_size;j++) {
      for (uint32_t k=0;k<f_size;k++) {
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
  for (uint32_t i=0;i<x_size;i++) {
    for (uint32_t j=0;j<y_size;j++) {
      for (uint32_t k=0;k<f_size;k++) {
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
  for (uint32_t i=0;i<n;i++) {
    output[i] = buffer[coords(x[i],y[i],f[i])];
  }
}

void skimage::set(uint32_t x, uint32_t y, uint32_t f, double val) {
  buffer[coords(x,y,f)] = val;
}
void skimage::set(uint32_t *x, uint32_t *y, uint32_t *f, double *val, uint32_t n) {
  for (uint32_t i=0;i<n;i++) {
    buffer[coords(x[i],y[i],f[i])] = val[i];
  }
}

double skimage::comp(uint32_t x, uint32_t y, uint32_t f, double val) {
  return val-buffer[coords(x,y,f)];
}
void skimage::comp(uint32_t *x, uint32_t *y, uint32_t *f, double *val, double *output, uint32_t n) {
  for (uint32_t i=0;i<n;i++) {
    output[i] = val[i] - buffer[coords(x[i],y[i],f[i])];
  }
}

double skimage::mean() {
  double runningtotal = 0.0;
  for (uint32_t i=0;i<totalsize;i++) {
    runningtotal += buffer[i];
  }
  return runningtotal/(double)totalsize;
}

double skimage::median() {
  vector <double> data;

  for (uint32_t i=0;i<totalsize;i++) {
    data.push_back(buffer[i]);
  }

  sort(data.begin(),data.end());

  return data[totalsize/2];
}

double skimage::max() {
  double max = 0.0;
  for (uint32_t i=0;i<totalsize;i++) {
    if (buffer[i]>max) { max=buffer[i]; }
  }
  return max;
}

double skimage::min() {
  double min = max();
  for (uint32_t i=0;i<totalsize;i++) {
    if (buffer[i]<min) { min=buffer[i]; }
  }
  return min;
}

void skimage::peak(uint32_t *x, uint32_t *y) {
  double max = 0.0;
  for (uint32_t v=0;v<y_size;v++) {
    for (uint32_t u=0;u<x_size;u++) {
      if (buffer[coords(u,v,0)]>max) {
        max=buffer[coords(u,v,0)];
        *x = u;
        *y = v;
      }
    }
  }
}

void skimage::unnan() {
  for (uint32_t i=0;i<totalsize;i++) {
    if (buffer[i]!=buffer[i]) {
      buffer[i] = 0;
    }
  }
}

double skimage::stdev(double mu) {
  double runningtotal = 0.0;
  for (uint32_t i=0;i<totalsize;i++) {
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
  for (uint32_t i=0;i<totalsize;i++) {
    dev = buffer[i]-mu0;
    if ((dev*dev)<cutoff) {
      mu += buffer[i]*pbeam.get(i);
      n+=pbeam.get(i);
    }
  }
  mu /= n;

  sig = 0;
  for (uint32_t i=0;i<totalsize;i++) {
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
  for (uint32_t i=0;i<totalsize;i++) {
    dev = buffer[i]-mu0;
    if ((dev*dev)<cutoff) {
      mu += buffer[i];
      n++;
    }
  }
  mu /= (double)n;

  sig = 0;
  for (uint32_t i=0;i<totalsize;i++) {
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
  for (uint32_t y=0;y<y_size;y++) {
    for (uint32_t x=0;x<x_size;x++) {
      logfile << " " << buffer[coords(x,y,0)];
    }
    logfile << endl;
  }

  for (uint32_t y=y_size-psize;y>0;y-=psize) {
    for (uint32_t x=0;x<x_size;x+=psize) {
      double pixel = 0.0;
      for (uint32_t u=0;u<psize;u++) {
        for (uint32_t v=0;v<psize;v++) {
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

  for (uint32_t f=0;f<f_size;f++) {
    for (uint32_t u=0;u<nx;u++) {
      for (uint32_t v=0;v<ny;v++) {
        uint32_t position = f+u*f_size+v*nx*f_size;
        uint32_t source = f+(u+dx)*f_size+(v+dy)*y_size*f_size;
        newbuffer[position] = oldbuffer[source];
      }
    }
  }

  x_0 += x_delt*(double)dx;
  y_0 += y_delt*(double)dy;

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

  for (uint32_t f=0;f<f_size;f++) {
    for (uint32_t u=0;u<nx;u++) {
      for (uint32_t v=0;v<ny;v++) {
        uint32_t position = f+u*f_size+v*nx*f_size;
        if (u<dx) { newbuffer[position]=contents; continue; }
        if (u>(nx-dx-1)) { newbuffer[position]=contents; continue; }
        if (v<dy) { newbuffer[position]=contents; continue; }
        if (v>(ny-dy-1)) { newbuffer[position]=contents; continue; }
        newbuffer[position] = oldbuffer[coords(u-dx,v-dy,f)];
      }
    }
  }

  x_0 -= x_delt*(double)dx;
  y_0 -= y_delt*(double)dy;

  buffer = newbuffer;
  x_size = nx;
  y_size = ny;

  totalsize = x_size*y_size*f_size;
  fxs = x_size*f_size;

  delete oldbuffer;
}

double *makekern(uint32_t ax, uint32_t ay, double pa, double maj, double min) {
  double *result;
  uint32_t cx = ax*0.5;
  uint32_t cy = ay*0.5;
  double rsq, flux, u, v;

  maj *= maj;
  min *= min;

  result = new double(ax*ay); 
  for (uint32_t y=0;y<ay;y++) {
    for (uint32_t x=0;x<ax;x++) {
      u = (x-cx)*cos(pa) - (y-cy)*sin(pa);
      v = (x-cx)*sin(pa) + (y-cy)*cos(pa);
      rsq = u*u/maj + v*v/min;
      if (rsq>1) {
        flux = -1;
      } else { 
        flux = exp(-rsq/2.);
      }
      result[y*ax+x] = flux;
    }
  }  

  return result;
}

double skimage::deconobject(skimage &other, skimage &pbeam, uint32_t *rax, uint32_t *ray, double *pa, double *min, double *maj, double *flux, uint64_t n) {
  double result = 0;
  vector <double> ourflux,ax,ay;
  uint32_t cx = other.getxsize()*0.5;
  uint32_t cy = other.getysize()*0.5;
  uint32_t kx, ky;
  double *kern;

  for (uint64_t i=0;i<n;i++) {
    //if (maj[n]>10) { maj[n]=10; }
    //if (min[n]>10) { min[n]=10; }
    kx = sqrt(-2.*maj[n]*maj[n]*log(noiseLevel/flux[n]));
    ky = kx;
    kern = makekern(2*kx, 2*ky, pa[n], maj[n], min[n]);
    for (uint32_t y=0;y<ky*2;y++) {
      for (uint32_t x=0;x<kx*2;x++) {
        if (kern[y*2*kx+x]>-1) {
          ourflux.push_back(flux[n]*kern[y*2*kx+x]);
          ax.push_back(x+rax[n]-kx);
          ay.push_back(y+ray[n]-ky);
        }
      }
    delete kern;
    }
  }
  n = ourflux.size();

  for (uint64_t i=0;i<n;i++) {
    ourflux[i]*=pbeam.get(ax[i],ay[i],0);
  }

  for (uint64_t i=0;i<n;i++) {
    result -= 2.0*ourflux[i]*buffer[coords(ax[i],ay[i],0)];
    for (uint64_t k=0;k<n;k++) {
      uint32_t u = cx - (ax[k]-ax[i]);
      uint32_t v = cy - (ay[k]-ay[i]);
      result += ourflux[i]*ourflux[k]*other.get(u,v,0);
    }
  }

  return -0.5*result/(noiseLevel*noiseLevel);
}

//Basic evaluator
double skimage::deconv(skimage &other, twov *points, double *flux, uint64_t n) {
  double result = 0;

  uint32_t cx = other.getxsize()*0.5;
  uint32_t cy = other.getysize()*0.5;
  for (uint64_t i=0;i<n;i++) {
    result += 2*flux[i]*buffer[coords(points[i].x, points[i].y,0)];
    for (uint64_t j=0;j<n;j++) {
      uint32_t u = cx + points[j].x-points[i].x;
      uint32_t v = cy + points[j].y-points[i].y;
      result -= flux[i]*flux[j]*other.get(u,v,0);
    }
  }

  return result * 0.5/(noiseLevel*noiseLevel);
}

//Evalutator including primary beam
double skimage::deconv(skimage &other, skimage &pbeam, uint32_t *ax, uint32_t *ay, double *flux, uint64_t n) {
  double result = 0;
  double ourflux[n];
  uint32_t cx = other.getxsize()*0.5;
  uint32_t cy = other.getysize()*0.5;

  for (uint64_t i=0;i<n;i++) {
    ourflux[i] = flux[i]*pbeam.get(ax[i],ay[i],0);
  }

  for (uint64_t i=0;i<n;i++) {
    result -= 2.0*ourflux[i]*buffer[coords(ax[i],ay[i],0)];
    for (uint64_t k=0;k<n;k++) {
      uint32_t u = cx - (ax[k]-ax[i]);
      uint32_t v = cy - (ay[k]-ay[i]);
      result += ourflux[i]*ourflux[k]*other.get(u,v,0);
    }
  }

  return -0.5*result/(noiseLevel*noiseLevel);
}

//Evalutator using frequency information and primary beam
double skimage::deconv(skimage &other, skimage &pbeam, twov *points, double *flux, double *fmu, double *fsig, uint64_t n) {
  double result = 0;

  uint32_t cx = other.getxsize()*0.5;
  uint32_t cy = other.getysize()*0.5;
  uint32_t cf = other.getfsize();

  for (uint32_t k=0;k<cf;k++) {
    for (uint64_t i=0;i<n;i++) {
      double fi = flux[i]*exp(-(k-fmu[i])*(k-fmu[i])/(fsig[i]*fsig[i]))*pbeam.get(points[i].x,points[i].y,k);
      result += 2*fi*buffer[coords(points[i].x, points[i].y,k)];
      for (uint64_t j=0;j<n;j++) {
        double fj = flux[j]*exp(-(k-fmu[j])*(k-fmu[j])/(fsig[j]*fsig[j]))*pbeam.get(points[j].x,points[j].y,k);
        uint32_t u = cx + points[j].x-points[i].x;
        uint32_t v = cy + points[j].y-points[i].y;
        result -= fi*fj*other.get(u,v,k);
      }
    }
  }

  return result * 0.5/(noiseLevel*noiseLevel);
}

void skimage::subbeam(skimage &other, twov points, double flux) {
  uint32_t cx = other.getxsize()*0.5;
  uint32_t cy = other.getysize()*0.5;

  for (uint32_t y=0;y<y_size;y++) {
    for (uint32_t x=0;x<x_size;x++) {
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

  for (uint32_t y=0;y<y_size;y++) {
    for (uint32_t x=0;x<x_size;x++) {
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

void skimage::setgrid(uint8_t axis, uint32_t crpix, double crval, double cdelt) {
  double zeroval = crval-((double)crpix*cdelt);
  if (axis==0) {
    x_0 = zeroval;
    x_delt = cdelt;
  } else {
    y_0 = zeroval;
    y_delt = cdelt;
  }
}

double skimage::real(uint8_t axis, double parval) {
  if (axis==0) {
    return parval*x_delt + x_0;
  }
  return parval*y_delt + y_0;
}
