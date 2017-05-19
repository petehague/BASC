#!/usr/bin/env python

import pymc
import model

S = pymc.MCMC(model, db='ram')
S.sample(iter=1000,burn=500)
#pymc.Matplot.plot(S)

output = open("trace.txt", "w")
output.write("x y flux\n")
for i in range(0,model.natoms):
    x = S.trace("xpos{}".format(i))[:]
    y = S.trace("ypos{}".format(i))[:]
    f = S.trace("flux{}".format(i))[:]
    for j in range(0,len(x)):
        output.write("{} {} {}\n".format(x[j],y[j],f[j]))
output.close()
