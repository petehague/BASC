#include <iostream>
#include <vector>
#include <cinttypes>
#include <fstream>
#include <sstream>
#include <chrono>

#include "../include/skimage.hpp"
#include "../bayesys/bayesys3.h"

skimage dirtyMap, dirtyBeam, primaryBeam;
skimage flatDirtyMap, flatDirtyBeam, flatPrimaryBeam;
double sigmasq, fluxscale, freqscale;

uint32_t imagesize, imagedepth;

bool cubeSwitch = false;

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
};

void setup(string mapfile, string psffile, string pbcorfile) {
  fstream imagefile;
  uint16_t y=0;
  uint32_t beamsize;

  dirtyMap.init(imagesize,imagesize,imagedepth);
  primaryBeam.init(imagesize,imagesize,imagedepth);

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

  dirtyBeam.init(beamsize,beamsize,imagedepth);

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


  sigmasq = dirtyMap.noise(primaryBeam);

  if (beamsize <= imagesize) {
    // cout << "Beam size " << beamsize << " padded out to " << beamsize+imagesize << endl;
    // dirtyBeam.pad(imagesize,imagesize,0);
    primaryBeam.crop(beamsize/2, beamsize/2);
    dirtyMap.crop(beamsize/2, beamsize/2);
    cout << "Image size " << imagesize << " cropped to " << beamsize/2 << endl;
    imagesize = beamsize/2;
  }

  //dirtyBeam.scan(32,dirtyBeam.max()/1000,dirtyBeam.max()/10000);

  //sigmasq = dirtyMap.noise(primaryBeam);
  dirtyMap.setnoise(sigmasq);
  cout << "Noise = " << sigmasq << " Jy/Beam" << endl;
  sigmasq *= sigmasq;

  //fluxscale = dirtyMap.fluxscale(dirtyBeam);
  fluxscale = sqrt(sigmasq);
  cout << "Flux Scale = " << fluxscale << endl;
  //fluxscale = 1;

  freqscale = 1./(double)imagedepth;
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
      cout << Common->cool << " " << UC->burnin << endl;
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
    }
  }

  if (UC->nmodels>UC->maxmodels) {
    return 1;
  } else {
    return 0;
  }
}
}

int main(int argc, char **argv) {
  fstream logfile, chainfile;
  chrono::system_clock::time_point startpoint;
  double elapsed, noise, median, atoms, datoms;
  uint32_t cx, cy, code;
  CommonStr Common;
  ObjectStr *Members;
  UserCommonStr UserCommon[1];

  if (argc<4) {
    cerr << "Too few arguments" << endl;
    return 0;
  }

  imagesize = stoi(argv[4]);
  cout << "Map size: " << imagesize << endl;
  imagedepth = stoi(argv[5]);
  cout << "Number of frequency channels: " << imagedepth << endl;
  if (imagedepth>1) { cubeSwitch = true; }

  setup(argv[1],argv[2],argv[3]);

  noise = dirtyMap.noise();
  median = dirtyMap.median();

  if (imagedepth>1) {
    Common.Ndim = 5;
  } else {
    Common.Ndim = 3;
  }

  Common.MinAtoms = 1;
  Common.MaxAtoms = 0;
  Common.Alpha = 1;
  Common.Valency = 0;
  Common.Iseed = 4321;
  Common.ENSEMBLE = 10;
  Common.Method = -1;
  Common.Rate = 1;
  Common.UserCommon = (void *)UserCommon;
  UserCommon->nmodels = 0;
  UserCommon->maxmodels = 300;
  UserCommon->burnin = 0;
  UserCommon->cool = 0;

  Members = new ObjectStr[Common.ENSEMBLE];

  startpoint = chrono::system_clock::now();
  dirtyMap.setnoise(9e-5);
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
    chainfile << UserCommon->x[i] << " " << UserCommon->y[i] << " " << UserCommon->F[i];
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
