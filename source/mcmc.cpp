#include <iostream>
#include <vector>
#include <cinttypes>
#include <fstream>
#include <sstream>

#include "../include/skimage.hpp"
#include "../bayesys/bayesys3.h"

skimage dirtyMap, dirtyBeam, primaryBeam;
double sigmasq, fluxscale;

uint32_t imagesize;

struct UserCommonStr {
  uint32_t nmodels;
  uint32_t maxmodels;
  vector <uint16_t> natoms;
  vector <double> x;
  vector <double> y;
  vector <double> F;
};

void setup(string mapfile, string psffile, string pbcorfile) {
  fstream imagefile;
  uint16_t y=0;

  dirtyMap.init(imagesize,imagesize,1);
  dirtyBeam.init(imagesize,imagesize,1);
  primaryBeam.init(imagesize,imagesize,1);

  imagefile.open(mapfile, fstream::in);
  for (string line; getline(imagefile, line); y++) {
    stringstream linestream(line);
    for (uint16_t x=0;x<imagesize;x++) {
      double value;
      linestream >> value;
      dirtyMap.set(y,x,0,value);
    }
  }
  imagefile.close();

  y=0;
  imagefile.open(psffile, fstream::in);
  for (string line; getline(imagefile, line); y++) {
    stringstream linestream(line);
    for (auto x=0;x<imagesize;x++) {
      double value;
      linestream >> value;
      dirtyBeam.set(y,x,0,value);
    }
  }
  imagefile.close();
  dirtyBeam.pad(imagesize,imagesize,0);

  y=0;
  imagefile.open(pbcorfile, fstream::in);
  for (string line; getline(imagefile, line); y++) {
    stringstream linestream(line);
    for (auto x=0;x<imagesize;x++) {
      double value;
      linestream >> value;
      primaryBeam.set(y,x,0,value);
    }
  }
  imagefile.close();

  primaryBeam.unnan();

  //primaryBeam.scan(16,primaryBeam.max()/2,primaryBeam.max()/4);

  sigmasq = dirtyMap.noise(primaryBeam);
  //sigmasq = 0.000086;
  dirtyMap.setnoise(sigmasq);
  cout << "Noise = " << sigmasq << " Jy/Beam" << endl;
  //sigmasq = 0.000083;
  sigmasq *= sigmasq;

  fluxscale = dirtyMap.fluxscale(dirtyBeam);
  fluxscale = sqrt(sigmasq);
  cout << "Flux Scale = " << fluxscale << endl;
  //fluxscale = 1;
}

extern "C" {
int UserBuild(double *like, CommonStr *Common, ObjectStr* Member, int natoms, int dummy) {
  vector <twov> points;
  vector <double> flux;
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
  }

  *like = dirtyMap.deconv(dirtyBeam,primaryBeam,&points[0],&flux[0],flux.size());

  //*like = dirtyMap.deconv(dirtyBeam,&points[0],&flux[0],flux.size());

  return 1;
}

int UserMonitor(CommonStr *Common, ObjectStr *Members) {
  double **Cube;// = Members->Cubes;
  UserCommonStr *UC = (UserCommonStr *)Common->UserCommon;

  if (Common->cool > 1.) Common->cool = 1.;
  if (Common->cool < 1.) return 0;

  UC->natoms.push_back(Members[0].Natoms);
  if (UC->nmodels == 0) {
    cout << "Burn in complete" << endl;
  }
  UC->nmodels += 1;
  for (auto k=0;k<Common->ENSEMBLE;k++) {
    Cube = Members[k].Cubes;
    for (auto i=0;i<Members[k].Natoms;i++) {
      UC->x.push_back(Cube[i][0]);
      UC->y.push_back(Cube[i][1]);
      UC->F.push_back(Cube[i][2]);
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

  setup(argv[1],argv[2],argv[3]);

  noise = dirtyMap.noise();
  median = dirtyMap.median();

  Common.Ndim = 3;
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

  Members = new ObjectStr[Common.ENSEMBLE];

  startpoint = chrono::system_clock::now();
  dirtyMap.setnoise(9e-5);
  code = BayeSys3(&Common,Members);
  elapsed = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(chrono::system_clock::now()-startpoint).count() )/1000.0;
  cout << "Chain time: " << elapsed << "s" << endl;

  cout << "Return code: " << code << endl;
  cout << "Evidence: " << Common.Evidence << endl;

  chainfile.open("chain.txt", fstream::out);
  for (auto i=0;i<UserCommon->x.size();i++) {
    chainfile << UserCommon->x[i] << " " << UserCommon->y[i] << " " << UserCommon->F[i] << endl;
    dirtyMap.subtract(twov(UserCommon->x[i],UserCommon->y[i]), UserCommon->F[i]/UserCommon->x.size());
  }
  chainfile.close();

  //Quick, dirty and wrong residual calculations
  //cout << "RMS Residuals: " << dirtyMap.residual(0,1) << " with noise " << dirtyMap.getnoise() << endl;

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
