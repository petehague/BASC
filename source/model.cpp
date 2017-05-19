/*
    Likelihood Template
    Author: Peter Hague
    Created: 22/07/14
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <cinttypes>
#include <vector>
#include <cmath>

#include "../include/skimage.hpp"
#include "../RainfallMCMC/include/agent.hpp"

using namespace std;

skimage dirtyMap, dirtyBeam;
double sigmasq, fluxscale;

class likelihood : public agent {
public:
    likelihood () {
        std::cout << "Created skimage model" << std::endl;
        dirtyBeam.init(512,512,1);
        dirtyMap.init(512,512,1);
    }

    void setup(options *o) {
        fstream imagefile;
        uint16_t y=0;

        imagefile.open("../NGC1808_cont.dirty.txt", fstream::in);

        for (string line; getline(imagefile, line); y++) {
          stringstream linestream(line);
          for (uint16_t x=0;x<512;x++) {
            double value;
            linestream >> value;
            dirtyMap.set(y,x,0,value);
          }
        }
        imagefile.close();

        y=0;
        imagefile.open("../NGC1808_cont.dirty.psf.txt", fstream::in);
        for (string line; getline(imagefile, line); y++) {
          stringstream linestream(line);
          for (auto x=0;x<512;x++) {
            double value;
            linestream >> value;
            dirtyBeam.set(y,x,0,value);
          }
        }
        dirtyBeam.pad(512,512,0);

        dirtyMap.scan(16,dirtyMap.mean(),dirtyMap.noise());
        dirtyBeam.scan(32,dirtyBeam.mean(),dirtyBeam.noise());

        sigmasq = dirtyMap.noise();
        cout << "Noise = " << sigmasq << " Jy/Beam" << endl;
        //sigmasq = 0.000083;
        sigmasq *= sigmasq;
        sigmasq = 1e-4;
        //sigmasq *= 1e4;

        fluxscale = dirtyMap.fluxscale(dirtyBeam);
        cout << "Flux Scale = " << fluxscale << endl;
        fluxscale = 1;

    }

    double eval(double *model) {
        vector<double> flux;
        vector<twov> points;
        double result;

        points.push_back(twov(model[0],model[1]));
        //flux.push_back(fluxscale*exp(model[2]));
        //return (model[0]-256)*(model[0]-256) + (model[1]-256)*(model[1]-256);
        flux.push_back(1);
        result = (0.5/sigmasq) * dirtyMap.deconv(dirtyBeam,&points[0],&flux[0],1);
        return result;
    }
};

REGISTERAGENT(likelihood)
