#include <iostream>
#include <vector>
#include <cinttypes>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>

#include "../include/skimage.hpp"
#include "../bayesys/bayesys3.h"
#include "../include/options.hpp"

skimage dirtyMap, dirtyBeam, primaryBeam;
skimage flatDirtyMap, flatDirtyBeam, flatPrimaryBeam;
skimage altMap, altBeam;
double sigma, fluxscale, freqscale;

uint32_t imagesize, imagedepth, modelIndex;

optDict options;

bool cubeSwitch = false;
bool altSwitch = false;

vector <uint32_t> sourcex;
vector <uint32_t> sourcey;
vector <double> sourceflux;
double sourcewidth;
double errscale;

struct UserCommonStr {
  uint32_t nmodels;
  uint32_t maxmodels;
  uint32_t burnin;
  double cool;
  vector <uint16_t> natoms;
  vector <double> x;
  vector <double> y;
  vector <double> F;
  vector <double> fmu;
  vector <double> fsig;
  vector <uint32_t> modelIndex;
};

void addSource(skimage map, uint32_t x, uint32_t y, double f, double sig) {
  uint32_t mx = map.getxsize();
  uint32_t my = map.getysize();

  //TODO Normalise flux

  for (auto v=0;v<my;v++) {
    for (auto u=0;u<mx;u++) {
      double rsq = (u-x)*(u-x)+(v-y)*(v-y);
      double flux = f*exp(-(rsq/(sig*sig)));
      map.add(twov(u,v), flux);
    }
  }
}

void setup(string mapfile, string psffile, string pbcorfile, string metafile) {
  fstream imagefile;
  uint16_t y=0;
  uint32_t beamsize;
  double x0,dx,y0,dy;
  uint32_t xref, yref;

  modelIndex = 0;

  imagefile.open(metafile, fstream::in);
  imagefile >> x0 >> dx >> xref >> y0 >> dy >> yref;
  imagefile.close();

  dirtyMap.init(imagesize,imagesize,imagedepth);
  dirtyMap.setgrid(0,xref,x0,dx);
  dirtyMap.setgrid(1,yref,y0,dy);
  primaryBeam.init(imagesize,imagesize,imagedepth);
  primaryBeam.setgrid(0,xref,x0,dx);
  primaryBeam.setgrid(1,yref,y0,dy);

  if (imagedepth>1) {
    flatDirtyMap.init(imagesize, imagesize, 1);
    flatPrimaryBeam.init(imagesize, imagesize, 1);
  }

  imagefile.open(mapfile, fstream::in);
  for (string line; getline(imagefile, line); y++) {
    stringstream linestream(line);
    for (uint16_t x=0;x<imagesize;x++) {
      double total = 0;
      for (uint16_t f=0;f<imagedepth;f++) {
        double value;
        linestream >> value;
        total += value;
        dirtyMap.set(y,x,f,value);
      }
      if (imagedepth>1) flatDirtyMap.set(y,x,0,1);
    }
  }
  imagefile.close();

  y=0;
  imagefile.open(pbcorfile, fstream::in);
  for (string line; getline(imagefile, line); y++) {
    stringstream linestream(line);
    for (auto x=0;x<imagesize;x++) {
      double total = 0;
      for (uint16_t f=0;f<imagedepth;f++) {
        double value;
        linestream >> value;
        total += value;
        primaryBeam.set(y,x,f,value);
      }
      if (imagedepth>1) flatPrimaryBeam.set(y,x,0,1);
    }
  }
  imagefile.close();

  primaryBeam.unnan();

  imagefile.open(psffile, fstream::in);
  beamsize = 0;
  for (string line; getline(imagefile, line); beamsize++);
  imagefile.close();
  imagefile.open(psffile, fstream::in);

  // Assumes equal beam size at this point
  dirtyBeam.init(beamsize,beamsize,imagedepth);
  dirtyBeam.setgrid(0,xref,x0,dx);
  dirtyBeam.setgrid(1,yref,y0,dy);

  if (imagedepth>1) {
    flatDirtyBeam.init(beamsize, beamsize, 1);
  }

  y = 0;
  for (string line; getline(imagefile, line); y++) {
    stringstream linestream(line);
    for (auto x=0;x<beamsize;x++) {
      double total = 0;
      for (uint16_t f=0;f<imagedepth;f++) {
        double value;
        linestream >> value;
        total += value;
        dirtyBeam.set(y,x,f,value);
      }
      if (imagedepth>1) flatDirtyBeam.set(y,x,0,1);
    }
  }
  imagefile.close();


  sigma = dirtyMap.noise(primaryBeam);

  if (beamsize <= imagesize) {
    // cout << "Beam size " << beamsize << " padded out to " << beamsize+imagesize << endl;
    // dirtyBeam.pad(imagesize,imagesize,0);
    primaryBeam.crop(beamsize/2, beamsize/2);
    dirtyMap.crop(beamsize/2, beamsize/2);
    cout << "Image size " << imagesize << " cropped to " << beamsize/2 << endl;
    imagesize = beamsize/2;
  }

  //dirtyBeam.scan(32,dirtyBeam.max()/1000,dirtyBeam.max()/10000);

  //sigma = dirtyMap.noise(primaryBeam);
  dirtyMap.setnoise(sigma);
  cout << "Noise = " << sigma << " Jy/Beam" << endl;
  sigma *= sigma;

  //fluxscale = dirtyMap.fluxscale(dirtyBeam);
  fluxscale = sqrt(sigma);
  //fluxscale = 1e-5;
  cout << "Flux Scale = " << fluxscale << endl;
  //fluxscale = 1;

  freqscale = 1./(double)imagedepth;

  //TODO add in metadata to this skimage
  if (altSwitch) {
    altMap.init(errscale, errscale, 1);
    altBeam.init(errscale, errscale, 1);
    for (auto i=0;i<sourcex.size();i++) {
      addSource(altMap, sourcex[i], sourcey[i], sourceflux[i], sourcewidth);
    }
    addSource(altBeam, errscale/2, errscale/2, 1, sourcewidth);
  }
}

extern "C" {
int UserBuild(double *like, CommonStr *Common, ObjectStr* Member, int natoms, int dummy) {
  vector <twov> points;
  vector <double> flux;
  vector <double> fmu;
  vector <double> fsig;
  double **Cube = Member->Cubes;
  UserCommonStr *UC = (UserCommonStr *)Common->UserCommon;

  if (natoms==0) {
    *like = -1e6;
    return 0;
  }

  for (auto i=0;i<natoms;i++) {
    double x = Cube[i][0]*double(imagesize);
    double y = Cube[i][1]*double(imagesize);
    points.push_back(twov(x,y));
    flux.push_back(fluxscale*Cube[i][2]/(1.-Cube[i][2]));
    if (cubeSwitch) {
      fmu.push_back(Cube[i][3]*freqscale);
      fsig.push_back(Cube[i][4]*freqscale);
    }
  }

  if (cubeSwitch) {
    *like = dirtyMap.deconv(dirtyBeam,primaryBeam,&points[0],&flux[0],&fmu[0],&fsig[0],flux.size());
    return 1;
  }

  *like = dirtyMap.deconv(dirtyBeam,primaryBeam,&points[0],&flux[0],flux.size());

  if (altSwitch) {
    *like = *like + altMap.deconv(altBeam,&points[0],&flux[0],flux.size());
  }

  //*like = dirtyMap.deconv(dirtyBeam,&points[0],&flux[0],flux.size());
  return 1;
}

int UserMonitor(CommonStr *Common, ObjectStr *Members) {
  double **Cube;// = Members->Cubes;
  UserCommonStr *UC = (UserCommonStr *)Common->UserCommon;

  if (Common->cool > 1.) Common->cool = 1.;
  if (Common->cool < 1.) {
    UC->burnin += 1;
    if (Common->cool>UC->cool) {
      //cout << Common->cool << " " << UC->burnin << endl;
      UC->cool += 0.1;
    }
    return 0;
  }

  UC->natoms.push_back(Members[0].Natoms);
  if (UC->nmodels == 0) {
    cout << "Burn in complete" << endl;
  }
  UC->nmodels += 1;
  for (auto k=0;k<Common->ENSEMBLE;k++) {
    Cube = Members[k].Cubes;
    for (auto i=0;i<Members[k].Natoms;i++) {
      UC->x.push_back(Cube[i][0]*double(imagesize));
      UC->y.push_back(Cube[i][1]*double(imagesize));
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
    return 1;
  } else {
    return 0;
  }
}

}

int main(int argc, char **argv) {
  fstream logfile, chainfile,sourcefile;
  chrono::system_clock::time_point startpoint;
  double elapsed, noise, median, atoms, datoms;
  int32_t cx, cy, code;
  CommonStr Common;
  ObjectStr *Members;
  UserCommonStr UserCommon[1];

  if (argc<5) {
    cerr << "Too few arguments" << endl;
    return 0;
  }

  imagesize = stoi(argv[5]);
  cout << "Map size: " << imagesize << endl;
  imagedepth = stoi(argv[6]);
  cout << "Number of frequency channels: " << imagedepth << endl;
  if (imagedepth>1) { cubeSwitch = true; }

  setup(argv[1],argv[2],argv[3],argv[4]);

  noise = dirtyMap.noise();
  median = dirtyMap.median();

  if (imagedepth>1) {
    Common.Ndim = 5;
  } else {
    Common.Ndim = 3;
  }

  options.readFile("config.txt");
  //options.report();

  if (options.getint("Sources")>0) {
    sourcefile.open("sources.txt", fstream::in);
    for (string line; getline(sourcefile, line); ) {
        stringstream linereader(line);
        double x,y,f,s;
        linereader >> x >> y >> f >> s;
        sourcex.push_back(x);
        sourcey.push_back(y);
        sourceflux.push_back(f);
        sourcewidth = s;
    }
    altSwitch = true;
    sourcefile.close();
  }

  Common.MinAtoms = options.getint("MinAtoms");
  Common.MaxAtoms = options.getint("MaxAtoms");
  Common.Alpha = options.getint("Alpha");
  Common.Valency = options.getint("Valency");
  Common.Iseed =  options.getint("Iseed");
  Common.ENSEMBLE =  options.getint("Ensemble");
  Common.Method =  options.getint("Method");
  Common.Rate =  options.getint("Rate");
  Common.UserCommon = (void *)UserCommon;
  UserCommon->nmodels = 0;
  UserCommon->maxmodels =  options.getint("maxmodels");
  UserCommon->burnin = 0;
  UserCommon->cool = 0;

  Members = new ObjectStr[Common.ENSEMBLE];

  startpoint = chrono::system_clock::now();

  code = BayeSys3(&Common,Members);
  elapsed = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(chrono::system_clock::now()-startpoint).count() )/1000.0;
  cout << "Chain time: " << elapsed << "s" << endl;

  cout << "Number of models: " << UserCommon->nmodels << endl;
  logfile.open("info.txt", fstream::out);
  logfile << UserCommon->nmodels << endl;
  logfile.close();

  cout << "Return code: " << code << endl;
  cout << "Evidence: " << Common.Evidence << endl;

  chainfile.open("chain.txt", fstream::out);
  for (auto i=0;i<UserCommon->x.size();i++) {
    chainfile << fixed << setprecision(9) << UserCommon->x[i] << " " << UserCommon->y[i] << " " << UserCommon->F[i] << " " << UserCommon->modelIndex[i];
    if (cubeSwitch) {
      chainfile << " " << UserCommon->fmu[i] << " " << UserCommon->fsig[i];
    }
    chainfile << endl;
  }
  chainfile.close();

  chainfile.open("chainreal.txt", fstream::out);
  for (auto i=0;i<UserCommon->x.size();i++) {
    chainfile << fixed << setprecision(9) << dirtyMap.real(0,UserCommon->x[i]) << " " << dirtyMap.real(1,UserCommon->y[i]) << " " << UserCommon->F[i] << " " << UserCommon->modelIndex[i];
    if (cubeSwitch) {
      chainfile << " " << UserCommon->fmu[i] << " " << UserCommon->fsig[i];
    }
    chainfile << endl;
  }
  chainfile.close();

  logfile.open("log.txt", fstream::out);
  atoms = 0;
  for (auto a : UserCommon->natoms) {
    logfile << a << endl;
    atoms += a;
    datoms += a*a;
  }
  logfile.close();

  atoms /= UserCommon->natoms.size();
  datoms = sqrt( (datoms - atoms*atoms*UserCommon->natoms.size()) / (UserCommon->natoms.size()/1) );

  cout << "Atoms: " << atoms << " +- " << datoms << endl;
}
