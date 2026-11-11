#!/usr/bin/env python3

#############################################################################
#
# MODULE:       r.avaflow
#
# AUTHORS:      Martin Mergili and Shiva P. Pudasaini
# CONTRIBUTORS: Massimiliano Alvioli, Matthias Benedikt, Emmanuel Delage,
#               Wolfgang Fellin, Jan-Thomas Fischer, Sigridur S. Gylfadottir,
#               Andreas Huber, Ivan Marchesini, Markus Metz, Markus Neteler,
#               Alexander Ostermann, Matthias Rauter
#
# PURPOSE:      The mass flow simulation tool
#
# COPYRIGHT:    (c) 2013 - 2025 by the authors
#               (c) 2020 - 2025 by the University of Graz
#               (c) 2013 - 2021 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (c) 2014 - 2020 by the University of Bonn
#               (c) 2020 - 2025 by the Technical University of Munich
#               (c) 2000 - 2025 by the GRASS Development Team
#               (c) 1993 - 2025 by the R Development Core Team
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################

#%module
#% description: The mass flow simulation tool
#% keywords: Raster
#% keywords: Landslide
#% keywords: Numerical simulation
#%end

#%flag
#% key: a
#% description: Map plots of pressure and kinetic energy
#% guisection: flags
#%end

#%flag
#% key: e
#% description: Model execution
#% guisection: flags
#%end

#%flag
#% key: k
#% description: Keep result GRASS raster maps
#% guisection: flags
#%end

#%flag
#% key: m
#% description: Multiple model runs
#% guisection: flags
#%end

#%flag
#% key: t
#% description: Map plots of impact wave or tsunami height
#% guisection: flags
#%end

#%flag
#% key: v
#% description: Evaluation and visualization
#% guisection: flags
#%end

#%option
#% key: prefix
#% type: string
#% description: Prefix for output files and folders
#% required: no
#% multiple: no
#%end

#%option
#% key: cores
#% type: string
#% description: Number of cores to be used for multiple model runs
#% required: no
#% multiple: no
#%end

#%option
#% key: cellsize
#% type: string
#% description: Cell size in metres
#% required: no
#% multiple: no
#%end

#%option
#% key: phases
#% type: string
#% description: Number of phases (1 = single-phase model, 3 = multi-phase model)
#% required: no
#% multiple: no
#%end

#%option
#% key: gravity
#% type: string
#% description: Gravity (m/s2)
#% required: no
#% multiple: no
#%end

#%option
#% key: limiter
#% type: string
#% description: Numerical limiter (1 = Minmod limiter, 2 = Superbee limiter, 3 = Woodward limiter, 4 = van Leer limiter)
#% required: no
#% multiple: no
#%end

#%option
#% key: clayers
#% type: string
#% description: Layer control (0 = of, 1 = on)
#% required: no
#% multiple: no
#%end

#%option
#% key: ctopo
#% type: string
#% description: Input and output of flow heights and velocities (0=vertical, 1=slope-normal)
#% required: no
#% multiple: no
#%end

#%option
#% key: csurface
#% type: string
#% description: Surface control (0 = off, 1 = edge control, 2 = force balancing)
#% required: no
#% multiple: no
#%end

#%option
#% key: centrainment
#% type: string
#% description: // Entrainment control (0 = off, 1=on)
#% required: no
#% multiple: no
#%end

#%option
#% key: cstopping
#% type: string
#% description: // Stopping control (0 = off, 1-3 = type of approach)
#% required: no
#% multiple: no
#%end

#%option
#% key: cdispersion
#% type: string
#% description: // Dispersion control (0 = off, 1 = on)
#% required: no
#% multiple: no
#%end

#%option
#% key: cmelt
#% type: string
#% description: // Melt control (0 = off, 1 = on)
#% required: no
#% multiple: no
#%end

#%option
#% key: aoicoords
#% type: string
#% description: Set of coordinates delineating area of interest (N,S,E,W)
#% required: no
#% multiple: yes
#%end

#%option
#% key: elevation
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of elevation raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hrelease
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of release height raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: rhrelease1
#% type: string
#% description: Ratio(s) of PHASE 1 release height (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: vhrelease
#% type: string
#% description: Variation of release height (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: hrelease1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 1 release height raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hrelease2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 2 release height raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hrelease3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 3 release height raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: trelease
#% type: string
#% gisprompt: old,raster,dcell
#% description: Release time raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: trelstop
#% type: string
#% gisprompt: old,raster,dcell
#% description: Release hydrograph stop time raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: tslide
#% type: string
#% gisprompt: old,raster,dcell
#% description: Time of initial sliding raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: tstop
#% type: string
#% gisprompt: old,raster,dcell
#% description: Stopping time raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: zones
#% type: string
#% gisprompt: old,raster,cell
#% description: Zones raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: impactarea
#% type: string
#% gisprompt: old,raster,cell
#% description: Name of observed impact area raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hdeposit
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of observed deposit height raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: vinx
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of release velocity in x direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: viny
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of release velocity in y direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: vinx1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 1 release velocity in x direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: viny1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 1 release velocity in y direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: vinx2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 2 release velocity in x direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: viny2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 2 release velocity in y direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: vinx3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 3 release velocity in x direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: viny3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 3 release velocity in y direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hentrmax
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of maximum height of entrainment raster map (comma-separated)
#% required: no
#% multiple: no
#%end

#%option
#% key: rhentrmax1
#% type: string
#% description: Ratio(s) of maximum height of PHASE 1 entrainment (comma-separated)
#% required: no
#% multiple: no
#%end

#%option
#% key: vhentrmax
#% type: string
#% description: Variation of maximum height of entrainment (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: hentrmax1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of maximum height of PHASE 1 entrainment raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hentrmax2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of maximum height of PHASE 2 entrainment raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hentrmax3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of maximum height of PHASE 3 entrainment raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: phi1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of internal friction angle of PHASE 1 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: phi2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of internal friction angle of PHASE 2 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: phi3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of internal friction angle of PHASE 3 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: delta1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of basal friction angle of PHASE 1 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: delta2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of basal friction angle of PHASE 2 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: delta3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of basal friction angle of PHASE 3 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: addfri1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of additional friction parameter of phase 1 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: addfri2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of additional friction parameter of phase 2 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: addfri3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of additional friction parameter of phase 3 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: coh1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of cohesion of phase 1 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: coh2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of cohesion of phase 2 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: coh3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of cohesion of phase 3 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ny1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of viscosity of PHASE 1 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ny2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of viscosity of PHASE 2 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ny3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of viscosity of PHASE 3 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: cdeform
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of deformation coefficient of PHASE 1 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: zfrag
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of fragmentation zones raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ambdrag
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of ambient drag coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: centr
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of entrainment coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ctrans12
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 1-PHASE 2 transformation coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ctrans13
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 1-PHASE 3 transformation coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ctrans23
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 2-PHASE 3 transformation coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: density
#% type: string
#% description: Density parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: friction
#% type: string
#% description: Friction parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: cohesion
#% type: string
#% description: Cohesion parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: viscosity
#% type: string
#% description: Viscosity parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: deformation
#% type: string
#% description: Deformation parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: slidepar
#% type: string
#% description: Parameters for initial sliding (exponent for weighting, coefficient for deformation)
#% required: no
#% multiple: yes
#%end

#%option
#% key: shearing
#% type: string
#% description: Strength of frictional energy loss in shearing
#% required: no
#% multiple: no
#%end

#%option
#% key: fragmentation
#% type: string
#% description: Fragmentation parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: ambient
#% type: string
#% description: Ambient drag coefficient
#% required: no
#% multiple: no
#%end

#%option
#% key: entrainment
#% type: string
#% description: Entrainment and stopping parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: drag
#% type: string
#% description: Drag parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: virtualmass
#% type: string
#% description: Virtual mass parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: transformation
#% type: string
#% description: Phase transformation parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: melting
#% type: string
#% description: Temperature evolution and ice melting parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: hydrograph
#% type: string
#% description: Path(es) to hydrograph file(s)
#% required: no
#% multiple: yes
#%end

#%option
#% key: hydrocoords
#% type: string
#% description: Pairs of coordinates of hydrograph(s) (x1, y1, x2, y2, ...)
#% required: no
#% multiple: yes
#%end

#%option
#% key: adaptograph
#% type: string
#% description: Path(es) to adaptograph file(s)
#% required: no
#% multiple: no
#%end

#%option
#% key: frictiograph
#% type: string
#% description: Path(es) to frictiograph file(s)
#% required: no
#% multiple: no
#%end

#%option
#% key: transformograph
#% type: string
#% description: Path(es) to transformograph file(s)
#% required: no
#% multiple: no
#%end

#%option
#% key: sampling
#% type: string
#% description: Type of parameter sampling (positive number = random, 0 = controlled, negative number = OAT)
#% required: no
#% multiple: no
#%end

#%option
#% key: cfl
#% type: string
#% description: CFL criterion, initial time step length (s)
#% required: no
#% multiple: yes
#%end

#%option
#% key: time
#% type: string
#% description: Time interval (s), stop for writing output (s)
#% required: no
#% multiple: yes
#%end

#%option
#% key: slomo
#% type: string
#% description: Parameters for slow motion (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: thresholds
#% type: string
#% description: Lower thresholds for display of flow height, velocity, kinetic energy, and pressure (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: profile
#% type: string
#% description: Pairs of coordinates for profile plots (x1,y1,x2,y2,...)
#% required: no
#% multiple: yes
#%end

#%option
#% key: ctrlpoints
#% type: string
#% description: Pairs of coordinates for control points (x1,y1,x2,y2,...)
#% required: no
#% multiple: yes
#%end

#%option
#% key: reftime
#% type: string
#% description: Reference time of reach for each control point (t1,t2,...)
#% required: no
#% multiple: yes
#%end

#%option
#% key: phexagg
#% type: string
#% description: Exaggeration of flow height in profile plots
#% required: no
#% multiple: no
#%end

#%option
#% key: orthophoto
#% type: string
#% description: Path to orthophoto tiff image
#% required: no
#% multiple: no
#%end

#%option
#% key: visualization
#% type: string
#% description: Parameters for visualization
#% required: no
#% multiple: yes
#%end

import grass.script as grass
from grass.script import core as grasscore
from grass.pygrass.modules import Module, MultiModule, ParallelModuleQueue
import math
from PIL import Image
import os
import queue
import random
import shutil
import subprocess
import sys
import threading
import time

def which(file):
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
            return path + "/" + file
    return None

for i in ["", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79"]:
    if not which("grass" + i) is None:
        bingrass = which("grass" + i)
        bingrass0 = "grass" + i

if not bingrass:
    grass.error("Please install GRASS 7 or GRASS 8.")
    grass.message(" ")
    sys.exit()

if not bingrass0 == "grass":
       
    queueLock = threading.Lock()
    workQueue = queue.Queue()

exitFlag = 0
ambvars = grass.gisenv()
locpath = ambvars.GISDBASE + "/" + ambvars.LOCATION_NAME
temppath = (
    ambvars.GISDBASE
    + "/"
    + ambvars.LOCATION_NAME
    + "/"
    + ambvars.MAPSET
    + "/.tmp/rtemp"
)

if bingrass0 == "grass":

    scriptpath = "$HOME/.grass8/addons/etc/r.avaflow.rcode"
    scriptpath2 = "$HOME/.grass8/addons/scripts"

else:

    scriptpath = "$HOME/.grass7/addons/etc/r.avaflow.rcode"
    scriptpath2 = "$HOME/.grass7/addons/scripts"

mainmapset = ambvars.MAPSET

def ErrorMessage(specify):
    grass.message(" ")
    grass.error("Please revise the " + specify + ".")
    grass.message(" ")
    sys.exit()

class myThread(threading.Thread):
    def __init__(self, name, q):
        threading.Thread.__init__(self)
        self.name = name
        self.q = q

    def run(self):
        process_data(self.name, self.q)

def StartBatch(jid):
    print("Executing model run %s" % jid)
    execute = "bash " + temppath + "/tmp" + str(jid) + "/batch" + str(jid)
    os.system(execute + " < /dev/null > " + temppath + "/out" + str(jid))
    print("Model run %s completed." % jid)
    return

def process_data(threadName, q):
    while not exitFlag:
        queueLock.acquire()
        if not workQueue.empty():
            time.sleep(0.25)
            data = q.get()
            queueLock.release()
            StartBatch(data)
        else:
            queueLock.release()

mstring = [
    "_hflow",
    "_tflow",
    "_pflow",
    "_basechange",
    "_vflow",
    "_htsun",
    "_treach",
]

def corrasc(corrname):

    os.rename(corrname + ".asc", corrname + "o.asc")

    fpressureo = open(corrname + "o.asc", "r")
    fpressure = open(corrname + ".asc", "w")

    fdata = ""
    for i in range(0, 5):
        fdatah = fpressureo.readline()
        fdata = fdata + "\t".join(fdatah.split()) + "\n"

    fdatah = fpressureo.readline()
    fdatai = fdatah.split()
    if fdatai[0] == "NODATA_value":
        fdata = fdata + fdatah
        fdataj = fdatai[1]
    else:
        fdata = fdata + "NODATA_value\t-9999\n" + fdatah
        fdataj = "-9999"

    fdata = fdata + "".join(fpressureo.readlines())
    fdata = fdata.replace("corner", "center").replace(fdataj, "-9999")

    fpressure.write(fdata)

    fpressureo.close()
    fpressure.close()

    os.system("rm -f " + corrname + "o.asc")
    os.system("rm -f " + corrname + ".asc.aux.xml")

def writeparam(jid, aflag, eflag, kflag, mflag, tflag, vflag, pf, cores, cellsize, phases, gravity, limiter, layers, ctopo, csurface, centrainment, cstopping, cdispersion, cmelt, aoicoords, 
   elevation, hrelease, rhrelease1, rhrls, vhrelease, vhrlx, hrelease1, hrelease2, hrelease3, trelease, trelstop, stoptime, tslide, zones, impactarea, hdeposit, 
   vinx, viny, vinx1, viny1, vinx2, viny2, vinx3, viny3, hentrmax, rhentrmax1, rhems, vhentrmax, vhemx, hentrmax1, hentrmax2, hentrmax3, 
   phi1, phi2, phi3, delta1, delta2, delta3, addfri1, addfri2, addfri3, coh1, coh2, coh3, ny1, ny2, ny3, cdeform, zfrag, ambdrag, centr, ctrans12, ctrans13, ctrans23, ortho_c1, ortho_c2, ortho_c3, 
   ortho_c4, ortho_c5, ortho_c6, hydrograph, hydrocoords, adaptograph, frictiograph, transformograph, gt, sampling, cfl, times, slomo, cvisc, cflux, thresholds, profile, ctrlpoints, 
   reftime, phexagg, orthophoto, lmax, nmesh, tint, tstop, visualization, ascpath):

   if not mflag:
       p1file = open(temppath + "/param1.txt", "w")
   else:
       p1file = open(temppath + "/param" + str(jid) + ".txt", "w")

   print(mainmapset, file=p1file)
   
   if not mflag:
       print("0", file=p1file)
   else:
       print("1", file=p1file)

   print("1", file=p1file)
   print("%s_" % pf, file=p1file)
   print("", file=p1file)
   print("%s_results/%s_ascii/" % (pf, pf), file=p1file)
   print("%s_results/%s_files/" % (pf, pf), file=p1file)
   print("%s_results/%s_aimec/" % (pf, pf), file=p1file)
   print("%s_results/%s_vr/" % (pf, pf), file=p1file)
   print("%s_results/%s_plots/" % (pf, pf), file=p1file)
   
   print(phases, file=p1file)
   
   if aflag:
       print("1", file=p1file)
   else:
       print("0", file=p1file)
   if tflag:
       print("1", file=p1file)
   else:
       print("0", file=p1file)
   print(gravity, file=p1file)
   print(limiter, file=p1file)
   print(layers, file=p1file)

   print(ctopo, file=p1file)
   print(csurface, file=p1file)
   print(centrainment, file=p1file)
   print(cstopping, file=p1file)
   print(cdispersion, file=p1file)
   print(cmelt, file=p1file)

   if elevation:
       print(elevation, file=p1file)
   else:
       print("None", file=p1file)            
   if hrelease:
       print(hrelease, file=p1file)
   elif hrelease1:
       print(hrelease1, file=p1file)
   else:
       print("None", file=p1file)
   if hrelease2:
       print(hrelease2, file=p1file)
   else:
       print("None", file=p1file)
   if hrelease3:
       print(hrelease3, file=p1file)
   else:
       print("None", file=p1file)
   if rhrelease1:
       print(rhrls, file=p1file)
   else:
       print("1.00", file=p1file)
   if vhrelease:
       print(vhrlx, file=p1file)
   else:
       print("1.00", file=p1file)
   if trelease:
       print(trelease, file=p1file)
   else:
       print("None", file=p1file)                    
   if trelstop:
       print(trelstop, file=p1file)
   else:
       print("None", file=p1file)                    
   if stoptime:
       print(stoptime, file=p1file)
   else:
       print("None", file=p1file)                    
   if tslide:
       print(tslide, file=p1file)
   else:
       print("None", file=p1file)
   if zones:
       print(zones, file=p1file)
   else:
       print("None", file=p1file)                     
   if impactarea:
       print(impactarea, file=p1file)
   else:
       print("None", file=p1file)                    
   if hdeposit:
       print(hdeposit, file=p1file)
   else:
       print("None", file=p1file)
   if vinx1:
       print(vinx1, file=p1file)
   else:
       print("None", file=p1file)
   if vinx2:
       print(vinx2, file=p1file)
   else:
       print("None", file=p1file)
   if vinx3:
       print(vinx3, file=p1file)
   else:
       print("None", file=p1file)
   if viny1:
       print(viny1, file=p1file)
   else:
       print("None", file=p1file)
   if viny2:
       print(viny2, file=p1file)
   else:
       print("None", file=p1file)
   if viny3:
       print(viny3, file=p1file)
   else:
       print("None", file=p1file)               
   if hentrmax == 1:
       print(hentrmax, file=p1file)
   elif hentrmax1:
       print(hentrmax1, file=p1file)
   else:
       print("None", file=p1file)
   if hentrmax2:
       print(hentrmax2, file=p1file)
   else:
       print("None", file=p1file)
   if hentrmax3:
       print(hentrmax3, file=p1file)
   else:
       print("None", file=p1file)
   if rhentrmax1:
       print(rhems, file=p1file)
   else:
       print("1.00", file=p1file)
   if vhentrmax:
       print(vhemx, file=p1file)
   else:
       print("1.00", file=p1file)
   if phi1:
       print(phi1, file=p1file)
   else:
       print("None", file=p1file)                    
   if phi2:
       print(phi2, file=p1file)
   else:
       print("None", file=p1file)                    
   if phi3:
       print(phi3, file=p1file)
   else:
       print("None", file=p1file)                    
   if delta1:
       print(delta1, file=p1file)
   else:
       print("None", file=p1file)
   if delta2:
       print(delta2, file=p1file)
   else:
       print("None", file=p1file)                    
   if delta3:
       print(delta3, file=p1file)
   else:
       print("None", file=p1file)
   if addfri1:
       print(addfri1, file=p1file)
   else:
       print("None", file=p1file)
   if addfri2:
       print(addfri2, file=p1file)
   else:
       print("None", file=p1file)
   if addfri3:
       print(addfri3, file=p1file)
   else:
       print("None", file=p1file)
   if coh1:
       print(coh1, file=p1file)
   else:
       print("None", file=p1file)
   if coh2:
       print(coh2, file=p1file)
   else:
       print("None", file=p1file)
   if coh3:
       print(coh3, file=p1file)
   else:
       print("None", file=p1file)
   if ny1:
       print(ny1, file=p1file)
   else:
       print("None", file=p1file)                    
   if ny2:
       print(ny2, file=p1file)
   else:
       print("None", file=p1file)                    
   if ny3:
       print(ny3, file=p1file)
   else:
       print("None", file=p1file)                    
   if cdeform:
       print(cdeform, file=p1file)
   else:
       print("None", file=p1file)  
   if zfrag:
       print(zfrag, file=p1file)
   else:
       print("None", file=p1file) 
   if ambdrag:
       print(ambdrag, file=p1file)
   else:
       print("None", file=p1file)
   if centr:
       print(centr, file=p1file)
   else:
       print("None", file=p1file)
   if ctrans12:
       print(ctrans12, file=p1file)
   else:
       print("None", file=p1file)                    
   if ctrans13:
       print(ctrans13, file=p1file)
   else:
       print("None", file=p1file)                    
   if ctrans23:
       print(ctrans23, file=p1file)
   else:
       print("None", file=p1file)

   if orthophoto:
   
       try: grass.parse_command("g.rename", rast=(ortho_c4, ortho_c1))
       except: test=0
       try: grass.parse_command("g.rename", rast=(ortho_c5, ortho_c2))
       except: test=0
       try: grass.parse_command("g.rename", rast=(ortho_c6, ortho_c3))
       except: test=0
       
       grass.run_command("r.out.gdal", input=ortho_c1, output=ascpath + ortho_c1 + ".asc", format="AAIGrid", overwrite=True)
       corrasc(ascpath + ortho_c1)     

       grass.run_command("r.out.gdal", input=ortho_c2, output=ascpath + ortho_c2 + ".asc", format="AAIGrid", overwrite=True)
       corrasc(ascpath + ortho_c2)

       grass.run_command("r.out.gdal", input=ortho_c3, output=ascpath + ortho_c3 + ".asc", format="AAIGrid", overwrite=True)
       corrasc(ascpath + ortho_c3)
   
       print(ortho_c1, file=p1file)
       print(ortho_c2, file=p1file)
       print(ortho_c3, file=p1file)
             
   else:
       print("None", file=p1file)
       print("None", file=p1file)
       print("None", file=p1file)

   if hydrocoords:
       print(hydrocoords, file=p1file)
   else:
       print("None", file=p1file)
       
   if hydrograph:
       print(hydrograph, file=p1file)
   else:
       print("None", file=p1file)

   if adaptograph:
       print(adaptograph, file=p1file)
   else:
       print("None", file=p1file)

   if frictiograph:
       print(frictiograph, file=p1file)
   else:
       print("None", file=p1file)

   if transformograph:
       print(transformograph, file=p1file)
   else:
       print("None", file=p1file)

   print(lmax, file=p1file)
   for l in range(0, lmax):
       print(round(gt[l], 10), file=p1file)

   print(thresholds[0], file=p1file)
   print(thresholds[1], file=p1file)
   print(thresholds[2], file=p1file)
   print(thresholds[3], file=p1file)
   print(thresholds[4], file=p1file)
   print(tint, file=p1file)
   print(tstop, file=p1file)
   print(slomo, file=p1file)
   print(cvisc, file=p1file)
   print(cflux, file=p1file)
   print(cfl[0], file=p1file)
   print(cfl[1], file=p1file)

   if profile:
       print(str(profile), file=p1file)
   else:
       print("None", file=p1file)

   if ctrlpoints:
       print(str(ctrlpoints), file=p1file)
   else:
       print("None", file=p1file)

   for i in range(0, len(visualization)):
       print(visualization[i], file=p1file)

   p1file.close()

def main():

    aflag = flags["a"]
    eflag = flags["e"]
    kflag = flags["k"]
    mflag = flags["m"]
    tflag = flags["t"]
    vflag = flags["v"]
    pf = options["prefix"]
    cores = options["cores"]
    cellsize = options["cellsize"]
    phases = options["phases"]
    gravity = options["gravity"]
    limiter = options["limiter"]
    layers = options["clayers"]
    ctopo = options["ctopo"]
    csurface = options["csurface"]
    centrainment = options["centrainment"]
    cstopping = options["cstopping"]
    cdispersion = options["cdispersion"]
    cmelt = options["cmelt"]
    aoicoords = options["aoicoords"]
    elevation = options["elevation"]
    hrelease = options["hrelease"]
    rhrelease1 = options["rhrelease1"]
    vhrelease = options["vhrelease"]
    hrelease1 = options["hrelease1"]
    hrelease2 = options["hrelease2"]
    hrelease3 = options["hrelease3"]
    trelease = options["trelease"]
    trelstop = options["trelstop"]
    tslide = options["tslide"]
    stoptime = options["tstop"]
    zones = options["zones"]
    impactarea = options["impactarea"]
    hdeposit = options["hdeposit"]
    vinx = options["vinx"]
    viny = options["viny"]
    vinx1 = options["vinx1"]
    viny1 = options["viny1"]
    vinx2 = options["vinx2"]
    viny2 = options["viny2"]
    vinx3 = options["vinx3"]
    viny3 = options["viny3"]
    hentrmax = options["hentrmax"]
    rhentrmax1 = options["rhentrmax1"]
    vhentrmax = options["vhentrmax"]
    hentrmax1 = options["hentrmax1"]
    hentrmax2 = options["hentrmax2"]
    hentrmax3 = options["hentrmax3"]
    phi1 = options["phi1"]
    phi2 = options["phi2"]
    phi3 = options["phi3"]
    delta1 = options["delta1"]
    delta2 = options["delta2"]
    delta3 = options["delta3"]
    addfri1 = options["addfri1"]
    addfri2 = options["addfri2"]
    addfri3 = options["addfri3"]
    coh1 = options["coh1"]
    coh2 = options["coh2"]
    coh3 = options["coh3"]
    ny1 = options["ny1"]
    ny2 = options["ny2"]
    ny3 = options["ny3"]
    cdeform = options["cdeform"]
    zfrag = options["zfrag"]
    ambdrag = options["ambdrag"]
    centr = options["centr"]
    ctrans12 = options["ctrans12"]
    ctrans13 = options["ctrans13"]
    ctrans23 = options["ctrans23"]
    orthophoto = options["orthophoto"]
    hydrograph = options["hydrograph"]
    hydrocoords = options["hydrocoords"]
    density = options["density"]
    friction = options["friction"]
    cohesion = options["cohesion"]
    viscosity = options["viscosity"]
    deformation = options["deformation"]
    slidepar = options["slidepar"]
    shearing = options["shearing"]
    fragmentation = options["fragmentation"]
    ambient = options["ambient"]
    entrainment = options["entrainment"]
    drag = options["drag"]
    virtualmass = options["virtualmass"]
    transformation = options["transformation"]
    melting = options["melting"]
    adaptograph = options["adaptograph"]
    frictiograph = options["frictiograph"]
    transformograph = options["transformograph"]
    sampling = options["sampling"]
    cfl = options["cfl"]
    times = options["time"]
    slomo = options["slomo"]
    thresholds = options["thresholds"]
    profile = options["profile"]
    ctrlpoints = options["ctrlpoints"]
    reftime = options["reftime"]
    phexagg = options["phexagg"]
    visualization = options["visualization"]

    if not pf:
        pf = "avf"

    print()
    print()
    print(
        "#############################################################################"
    )
    print("#")
    print("# MODULE:       r.avaflow 4.0 G Revision 7")
    print("#")
    print("# AUTHORS:      Martin Mergili and Shiva P. Pudasaini")
    print("# CONTRIBUTORS: Massimiliano Alvioli, Matthias Benedikt, Emmanuel Delage,")
    print("#               Wolfgang Fellin, Jan-Thomas Fischer, Sigridur S. Gylfadottir,")
    print("#               Andreas Huber, Ivan Marchesini, Markus Metz, Markus Neteler,")
    print("#               Alexander Ostermann, Matthias Rauter")
    print("#")
    print("# PURPOSE:      The mass flow simulation tool")
    print("# VERSION:      20250902 (2 September 2025)")
    print("#")
    print("# COPYRIGHT:    (c) 2013 - 2025 by the authors")
    print("#               (c) 2020 - 2025 by the University of Graz")
    print("#               (c) 2013 - 2021 by the BOKU University, Vienna")
    print("#               (c) 2015 - 2020 by the University of Vienna")
    print("#               (c) 2014 - 2020 by the University of Bonn")
    print("#               (c) 2020 - 2025 by the Technical University of Munich")
    print("#               (c) 2000 - 2025 by the GRASS Development Team")
    print("#               (c) 1993 - 2025 by the R Development Core Team")
    print("#")
    print("#               This program is free software under the GNU General Public")
    print("#               License (>=v2). Read the file COPYING that comes with GRASS")
    print("#               for details.")
    print("#")
    print("# PREFIX:       %s" % pf)
    print("#")
    print(
        "#############################################################################"
    )
    print()
    print()

    if not eflag and not vflag:

        eflag = "1"
        vflag = "1"

    if mflag:
        aflag = "1"

    if profile:
        if not phexagg:
            phexagg = "1"
        try:
            phexagg = float(phexagg)
        except ValueError:
            ErrorMessage("profile exaggeration value")
        phexagg = str(phexagg)

    if eflag:

        for inrast in [ elevation, hrelease, hrelease1, hrelease2, hrelease3, trelease, trelstop, tslide, stoptime, zones, impactarea, hdeposit, 
            vinx, viny, vinx1, viny1, vinx2, viny2, vinx3, viny3, hentrmax, hentrmax1, hentrmax2, hentrmax3, phi1, phi2, phi3,
            delta1, delta2, delta3, addfri1, addfri2, addfri3, coh1, coh2, coh3, ny1, ny2, ny3, cdeform, zfrag, ambdrag, centr, ctrans12, ctrans13, ctrans23 ]:

            if inrast and (not grasscore.find_file(inrast, element="cell", mapset=".") 
                and not grasscore.find_file(inrast, element="fcell", mapset=".") and not grasscore.find_file(inrast, element="dcell", mapset=".")):
                ErrorMessage("GRASS raster maps")

        nmesh = "1"

        if not phases:
            phases = "1"
        elif not phases == "1" and not phases == "3":
            ErrorMessage("number of phases")

        if not gravity:
            gravity = "9.81"

        if not limiter:
            limiter = "1"
        elif not limiter == "1" and not limiter == "2" and not limiter == "3" and not limiter == "4":
            ErrorMessage("numerical limiter")

        if not layers:
            layers = "0"
        elif not layers == "0" and not layers == "1" and not layers == "2":
            ErrorMessage("layer control")

        if not ctopo:
            ctopo = "0"
        elif not ctopo == "0" and not ctopo == "1":
            ErrorMessage("input and output control")

        if not csurface:
            csurface = "0"
        elif not csurface == "0" and not csurface == "1" and not csurface == "2" and not csurface == "3":
            ErrorMessage("surface control")

        if not centrainment:
            centrainment = "0"
        elif centrainment not in ("0", "1", "9", "-1"):
            ErrorMessage("entrainment control")

        if not cstopping:
            cstopping = "0"
        elif not cstopping == "0" and not cstopping == "1" and not cstopping == "2" and not cstopping == "3":
            ErrorMessage("stopping control")

        if not cdispersion:
            cdispersion = "0"
        elif not cdispersion == "0" and not cdispersion == "1" and not cdispersion == "2" and not cdispersion == "3":
            ErrorMessage("dispersion control")

        if not cmelt:
            cmelt = "0"
        elif not cmelt == "0" and not cmelt == "1":
            ErrorMessage("melt control")

        if not elevation:
            ErrorMessage("elevation raster map")

        if phases == "1" and not hydrograph and not hrelease:
            ErrorMessage("definition of flow release")
        elif phases == "3" and not hydrograph and not hrelease1 and not hrelease2 and not hrelease3 and not (hrelease and rhrelease1):
            ErrorMessage("definition of flow release")

        if phases == "3" and (hentrmax and not rhentrmax1):
            ErrorMessage("definition of entrainable height")

        if hydrograph and not hydrocoords:
            ErrorMessage("input hydrograph coordinates")

        if not thresholds:
            thresholds = "0.1,10000,10000,0.1,0.000001"
        thresholds = list(map(str, thresholds.split(",")))
        if not len(thresholds) == 5:
            ErrorMessage("number of thresholds")

        for ethresholds in thresholds:
            try:
                ethresholds = float(ethresholds)
            except ValueError:
                ErrorMessage("threshold parameter values")

        if mflag:
            if not sampling:
                sampling = "100"
            try:
                sampling = int(sampling)
            except ValueError:
                ErrorMessage("value provided for sampling")
            sampling=str(sampling)

        if not mflag:

            if phases == "1":

                if not density:
                    density = "2700"
                if not friction:
                    friction = "35,20,0.0"
                if not cohesion:
                    cohesion = "0.0"
                if not viscosity:
                    viscosity = "-7.0"
                if not deformation:
                    deformation = "1.0"
                if not slidepar:
                    slidepar = "0.0,0.0"
                if not shearing:
                    shearing = "0.0"
                if not fragmentation:
                    fragmentation = "0.0,0.0"
                if not ambient:
                    ambient = "0.0"
                if not entrainment:
                    entrainment = "-7.0,0.0"
                lmax = 15

            else:

                if not density:
                    density = "2700,1850,1000"
                if not friction:
                    friction = "35,17.5,0,20,10,0,0.0,0.0,0.0"
                if not cohesion:
                    cohesion = "0.0,0.0,0.0"
                if not viscosity:
                    viscosity = "-7.0,-7.0,-7.0"
                if not deformation:
                    deformation = "1.0,1.0,1.0"
                if not slidepar:
                    slidepar = "0.0,0.0,0.0,0.0,0.0,0.0"
                if not shearing:
                    shearing = "0.0"
                if not fragmentation:
                    fragmentation = "0.0,0.0"
                if not ambient:
                    ambient = "0.0"
                if not entrainment:
                    entrainment = "-7.0,0.0"
                if not drag:
                    drag = "1,3,1,0.1,1,1"
                if not virtualmass:
                    virtualmass = "10,0.12,1"
                if not transformation:
                    transformation = "0.0,0.0,0.0"
                if not melting:
                    melting = "0,0,0,0.2,0.5"
                lmax = 50

        elif int(sampling) > 0:

            if phases == "1":

                if not density:
                    density = "2700,2700"
                if not friction:
                    friction = "35,35,20,20,0.0,0.0"
                if not cohesion:
                    cohesion = "0.0,0.0"
                if not viscosity:
                    viscosity = "-7.0,-7.0"
                if not deformation:
                    deformation = "1.0,1.0"
                if not slidepar:
                    slidepar = "0.0,0.0,0.0,0.0"
                if not shearing:
                    shearing = "0.0,0.0"
                if not fragmentation:
                    fragmentation = "0.0,0.0,0.0,0.0"
                if not ambient:
                    ambient = "0.0,0.0"
                if not entrainment:
                    entrainment = "-7.0,-7.0,0.0,0.0"
                lmax = 15

            else:

                if not density:
                    density = "2700,2700,1850,1850,1000,1000"
                if not friction:
                    friction = "35,35,17.5,17.5,0,0,0,20,20,10,10,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0"
                if not cohesion:
                    cohesion = "0.0,0.0,0.0,0.0,0.0,0.0"
                if not viscosity:
                    viscosity = "-7.0,-7.0,-7.0,-7.0,-7.0,-7.0"
                if not deformation:
                    deformation = "1.0,1.0,1.0,1.0,1.0,1.0"
                if not slidepar:
                    slidepar = "0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0"
                if not shearing:
                    shearing = "0.0,0.0"
                if not fragmentation:
                    fragmentation = "0.0,0.0,0.0,0.0"
                if not ambient:
                    ambient = "0.0,0.0"
                if not entrainment:
                    entrainment = "-7.0,-7.0,0.0,0.0"
                if not drag:
                    drag = "1,1,3,3,1,1,0.1,0.1,1,1,1,1"
                if not virtualmass:
                    virtualmass = "10,10,0.12,0.12,1,1"
                if not transformation:
                    transformation = "0.0,0.0,0.0,0.0,0.0,0.0"
                if not melting:
                    melting = "0,0,0,0,0,0,0.2,0.2,0.5,0.5"
                lmax = 50

        else:

            if phases == "1":

                if not density:
                    density = "2700,2700,0"
                if not friction:
                    friction = "35,35,0,20,20,0,0.0,0.0,0"
                if not cohesion:
                    cohesion = "0.0,0.0,0"
                if not viscosity:
                    viscosity = "-7.0,-7.0,0"
                if not deformation:
                    deformation = "1.0,1.0,0"
                if not slidepar:
                    slidepar = "0.0,0.0,0,0.0,0.0,0"                    
                if not shearing:
                    shearing = "0.0,0.0,0"
                if not fragmentation:
                    fragmentation = "0.0,0.0,0,0.0,0.0,0"
                if not ambient:
                    ambient = "0.0,0.0,0"
                if not entrainment:
                    entrainment = "-7.0,-7.0,0,0.0,0.0,0"
                lmax = 15

            else:

                if not density:
                    density = "2700,2700,0,1850,1850,0,1000,1000,0"
                if not friction:
                    friction = "35,35,0,17.5,17.5,0,0,0,0,0,20,20,0,10,10,0,0,0,0,0,0.0,0.0,0,0.0,0.0,0,0.0,0.0,0"
                if not cohesion:
                    cohesion = "0.0,0.0,0,0.0,0.0,0,0.0,0.0,0"
                if not viscosity:
                    viscosity = "-7.0,-7.0,0,-7.0,-7.0,0,-7.0,-7.0,0"
                if not deformation:
                    deformation = "1.0,1.0,0,1.0,1.0,0,1.0,1.0,0"
                if not slidepar:
                    slidepar = "0.0,0.0,0,0.0,0.0,0,0.0,0.0,0,0.0,0.0,0,0.0,0.0,0,0.0,0.0,0"
                if not shearing:
                    shearing = "0.0,0.0,0"
                if not fragmentation:
                    fragmentation = "0.0,0.0,0,0.0,0.0,0"
                if not ambient:
                    ambient = "0.0,0.0,0"
                if not entrainment:
                    entrainment = "-7.0,-7.0,0,0.0,0.0,0"
                if not drag:
                    drag = "1,1,0,3,3,0,1,1,0,0.1,0.1,0,1,1,0,1,1,0"
                if not virtualmass:
                    virtualmass = "10,10,0,0.12,0.12,0,1,1,0"
                if not transformation:
                    transformation = "0.0,0.0,0,0.0,0.0,0,0.0,0.0,0"
                if not melting:
                    melting = "0,0,0,0,0,0,0,0,0,0.2,0.2,0,0.5,0.5,0"
                lmax = 50

        if phases == "1":
            flowparam = ( density + "," + friction + "," + cohesion + "," + viscosity + "," + deformation + "," + slidepar + "," + shearing + "," + fragmentation + "," + ambient + "," 
                + entrainment )
        else:
            flowparam = ( density + "," + friction + "," + cohesion + "," + viscosity + "," + deformation + "," + slidepar + "," + shearing + "," + fragmentation + "," + ambient + "," 
                + entrainment + "," + drag + "," + virtualmass + "," + transformation + "," + melting )

        flowparam = list(map(str, flowparam.split(",")))

        if not mflag and not len(flowparam) == lmax:
            ErrorMessage("number of flow parameters")
        elif mflag and int(sampling) > 0 and not len(flowparam) == 2 * lmax:
            ErrorMessage("number of flow parameters")
        elif mflag and int(sampling) <= 0 and not len(flowparam) == 3 * lmax:
            ErrorMessage("number of flow parameters")

        for eflowparam in flowparam:
            try:
                eflowparam = float(eflowparam)
            except ValueError:
                ErrorMessage("parameter values")

        if not cfl:
            cfl = "0.40,0.001"
        cfl = list(map(str, cfl.split(",")))
        if not len(cfl) == 2:
            ErrorMessage("number of cfl parameters")

        for ecfl in cfl:
            try:
                ecfl = float(ecfl)
            except ValueError:
                ErrorMessage("cfl values")

        if not times:
            times = "10,300"
        times = list(map(str, times.split(",")))
        tint = times[0]
        tstop = times[1]

        if mflag:
            tint = tstop

        for etimes in times:
            try:
                etimes = float(etimes)
            except ValueError:
                ErrorMessage("times")

        if reftime:

            reftime = list(map(str, reftime.split(",")))
            
            ctrlpts = list(map(str, ctrlpoints.split(",")))
            if not len(reftime) == len(ctrlpts)/2:
                ErrorMessage("number of reference times of reach")

        if not slomo:
            slomo = "1.0,1.0,1.0"

        elif slomo == "min":
            slomo = "60,1.0,1.0"

        elif slomo == "h":
            slomo = "3600,1.0,1.0"
            
        elif slomo == "d":
            slomo = "86400,1.0,1.0"
            
        elif slomo == "w":
            slomo = "604800,1.0,1.0"

        elif slomo == "mon":
            slomo = "2592000,1.0,1.0"

        elif slomo == "yr":
            slomo = "31536000,1.0,1.0"

        elif slomo == "g":
            slomo = "-31536000,1.0,1.0"

        lslomo = list(map(str, slomo.split(",")))
        slomo = lslomo[0]
        cvisc = lslomo[1]
        cflux = lslomo[2]

        if mflag:
            if not cores:
                cores = "8"
            try:
                cores = int(cores)
            except ValueError:
                ErrorMessage("number of cores")
            cores=str(cores)

        if not visualization:
            visualization = "0,0.1,5.0,5.0,1,100,2,-11000,9000,100,0.60,0.25,0.15,0.2,1.0,C:/Program Files/ParaView 5.11.0/bin/pvpython.exe,C:/Program Files/R/R-4.2.3/bin/Rscript.exe,None"
        visualization = list(map(str, visualization.split(",")))
        if not len(visualization) == 18:
            ErrorMessage("number of visualization parameters")

    print("1. PREPARING ENVIRONMENT.")

    os.system("rm -rf " + temppath)
    os.system("mkdir " + temppath)

    os.environ["GRASS_VERBOSE"] = "-1"
    grass.run_command("g.gisenv", set="GRASS_VERBOSE=-1")

    ortho_c1 = pf + "_ortho.red"
    ortho_c2 = pf + "_ortho.green"
    ortho_c3 = pf + "_ortho.blue"
    ortho_c4 = pf + "_ortho.1"
    ortho_c5 = pf + "_ortho.2"
    ortho_c6 = pf + "_ortho.3"
    if orthophoto: grass.run_command("r.in.gdal", flags="o", overwrite=True, input=orthophoto, output=pf + "_ortho")

    ascpath = pf + "_results/" + pf + "_ascii/"
    filepath = pf + "_results/" + pf + "_files/"

    if eflag:

        if os.path.exists(pf + "_results"):
            os.system("rm -rf " + pf + "_results")
        os.mkdir(pf + "_results")
        os.mkdir(pf + "_results/" + pf + "_files")
        os.mkdir(pf + "_results/" + pf + "_plots")
        if mflag: 
            os.mkdir(pf + "_results/" + pf + "_aimec")
            os.mkdir(pf + "_results/" + pf + "_aimec/depth")
            os.mkdir(pf + "_results/" + pf + "_aimec/pressure")
        os.mkdir(pf + "_results/" + pf + "_vr")
        os.mkdir(pf + "_results/" + pf + "_vr/data")
        if not mflag:
            os.mkdir(pf + "_results/" + pf + "_plots/" + pf + "_maps_timesteps")
        if not mflag:
            os.mkdir(pf + "_results/" + pf + "_plots/" + pf + "_profiles_timesteps")
        os.mkdir(ascpath)

        grass.run_command("g.region", flags="d")
        if cellsize:
            grass.run_command("g.region", flags="a", res=cellsize)
        if aoicoords:
            aoicoords = list(map(str, aoicoords.split(",")))
            if not len(aoicoords) == 4:
                ErrorMessage("number of coordinates of the area of interest")
            grass.run_command("g.region", flags="a", n=aoicoords[0], s=aoicoords[1], w=aoicoords[2], e=aoicoords[3])

        grass.mapcalc('"_aoi"=pow("%s",0)' % elevation, overwrite=True, quiet=True)
        cellsize = grass.raster_info("_aoi")["nsres"]
        rnorth = grass.raster_info("_aoi")["north"]
        rsouth = grass.raster_info("_aoi")["south"]
        rwest = grass.raster_info("_aoi")["west"]
        reast = grass.raster_info("_aoi")["east"]

        grass.run_command("g.region", flags="d")
        grass.run_command("g.region", flags="a", n=rnorth, s=rsouth, w=rwest, e=reast)
        grass.run_command("g.region", flags="a", res=cellsize)

        if not mflag:

            jid = 0
            vhrl = "1.00"
            vhem = "1.00"
            gt = [float(flowparam[0])]
            for l in range(1, lmax):
                gt.append(float(flowparam[l]))
                
            if rhrelease1:
           
                rhrls = rhrelease1
               
            else:
           
                rhrls = "-9999"

            if rhentrmax1:
           
                rhems = rhentrmax1
               
            else:
           
                rhems = "-9999"

            writeparam(jid, aflag, eflag, kflag, mflag, tflag, vflag, pf, cores, cellsize, phases, gravity, limiter, layers, ctopo, csurface, centrainment, cstopping, cdispersion, cmelt, aoicoords, 
                elevation, hrelease, rhrelease1, rhrls, vhrelease, vhrl, hrelease1, hrelease2, hrelease3, trelease, trelstop, stoptime, tslide, zones, impactarea, hdeposit, 
                vinx, viny, vinx1, viny1, vinx2, viny2, vinx3, viny3, hentrmax, rhentrmax1, rhems, vhentrmax, vhem, hentrmax1, hentrmax2, hentrmax3, 
                phi1, phi2, phi3, delta1, delta2, delta3, addfri1, addfri2, addfri3, coh1, coh2, coh3, ny1, ny2, ny3, cdeform, zfrag, ambdrag, centr, ctrans12, ctrans13, ctrans23, 
                ortho_c1, ortho_c2, ortho_c3, ortho_c4, ortho_c5, ortho_c6, hydrograph, hydrocoords,adaptograph, frictiograph, transformograph, gt, sampling, cfl, times, slomo, 
                cvisc, cflux, thresholds, profile, ctrlpoints, reftime, phexagg, orthophoto, lmax, nmesh, tint, tstop, visualization, ascpath)

            print("2. ROUTING FLOW.")

            start = time.time()

            os.environ["XINT"] = "1"
            grass.run_command("r.avaflow.main")

            stop = time.time()
            comptime = stop - start

            timefile = open(filepath + pf + "_time.txt", "w")
            timefile.write(str(comptime))
            timefile.close()

        else:
        
            print("2. ROUTING FLOW.")
            if not bingrass0 == "grass":
                print("")

            if int(sampling) > 0:
                nruns = int(sampling)
            elif int(sampling) < 0:
                nruns = int(sampling) * -1
            ncores = int(cores)

            if phases == "3" and hrelease:
                rhrls = list(
                    map(str, rhrelease1.split(","))
                )
                rhrlsmin = rhrls[0]
                if len(rhrls) == 3:
                    rhrlsmax = rhrls[1]
                    rhrls3 = rhrls[2]
                elif len(rhrls) == 2:
                    rhrlsmax = rhrls[1]
                elif len(rhrls) == 1:
                    rhrlsmax = rhrls[0]
                else:
                    ErrorMessage()

            if (hrelease or hrelease1 or hrelease2 or hrelease3) and vhrelease:
                vhrl = list(
                    map(str, vhrelease.split(","))
                )
                vhrlmin = vhrl[0]
                if len(vhrl) == 3:
                    vhrlmax = vhrl[1]
                    vhrl3 = vhrl[2]
                elif len(vhrl) == 2:
                    vhrlmax = vhrl[1]
                elif len(vhrl) == 1:
                    vhrlmax = vhrl[0]
                else:
                    ErrorMessage()
            else:
                vhrlmin = "1"
                vhrlmax = "1"
                vhrl3 = "1"

            if phases == "3" and hentrmax:
                rhems = list(
                    map(str, rhentrmax1.split(","))
                )
                rhemsmin = rhems[0]
                if len(rhems) == 3:
                    rhemsmax = rhems[1]
                    rhems3 = rhems[2]
                elif len(rhems) == 2:
                    rhemsmax = rhems[1]
                elif len(rhems) == 1:
                    rhemsmax = rhems[0]
                else:
                    ErrorMessage()

            if (hentrmax or hentrmax1 or hentrmax2 or hentrmax3) and vhentrmax:
                vhem = list(map(str, vhentrmax.split(",")))
                vhemmin = vhem[0]
                if len(vhem) == 2:
                    vhemmax = vhem[1]
                    vhem3 = vhem[2]
                elif len(vhem) == 2:
                    vhemmax = vhem[1]
                elif len(vhem) == 1:
                    vhemmax = vhem[0]
                else:
                    ErrorMessage()
            else:
                vhemmin = "1"
                vhemmax = "1"
                vhem3 = "1"

            if sampling == "0":

                njids = 2
                nruns0 = []
                nruns = 1

                if float(vhrlmin) == float(vhrlmax):
                    nruns0.append(1)
                else:
                    nruns0.append(int(vhrl3))
                nruns *= nruns0[0]

                if float(vhemmin) == float(vhemmax):
                    nruns0.append(1)
                else:
                    nruns0.append(int(vhem3))
                nruns *= nruns0[1]

                if phases == "3" and hrelease and rhrelease1:

                    njids += 1
                    if float(rhrlsmin) == float(rhrlsmax):
                        nruns0.append(1)
                    else:
                        nruns0.append(int(rhrls3))
                    nruns *= nruns0[2]

                if phases == "3" and hentrmax and rhentrmax1:

                    njids += 1
                    if float(rhemsmin) == float(rhemsmax):
                        nruns0.append(1)
                    else:
                        nruns0.append(int(rhems3))
                    nruns *= nruns0[njids - 1]

                for l in range(0, lmax):

                    if float(flowparam[3 * l]) == float(flowparam[3 * l + 1]):
                        nruns0.append(1)
                    else:
                        nruns0.append(int(flowparam[3 * l + 2]))
                    nruns *= nruns0[l + njids]

                jids = [0]
                for l in range(1, lmax + njids):
                    jids.append(0)

            elif int(sampling) < 0:

                nruns00 = nruns
                nruns = 0
                nvar = []
                njids = 2

                if not float(vhrlmin) == float(vhrlmax):
                    nvar.append(1)
                    nruns += nruns00
                else:
                    nvar.append(0)

                if not float(vhemmin) == float(vhemmax):
                    nvar.append(1)
                    nruns += nruns00
                else:
                    nvar.append(0)

                if phases == "3" and hrelease and rhrelease1:

                    njids += 1
                    if not float(rhrlsmin) == float(rhrlsmax):
                        nvar.append(1)
                        nruns += nruns00
                    else:
                        nvar.append(0)

                if phases == "3" and hentrmax and rhentrmax1:

                    njids += 1
                    if not float(rhemsmin) == float(rhemsmax):
                        nvar.append(1)
                        nruns += nruns00
                    else:
                        nvar.append(0)

                for l in range(0, lmax):
                    if not float(flowparam[3 * l]) == float(flowparam[3 * l + 1]):
                        nvar.append(1)
                        nruns += nruns00
                    else:
                        nvar.append(0)

                lnrun = 0
                ltest = 0

            process_queue = ParallelModuleQueue(nprocs=ncores)

            for jid in range(1, nruns + 1):

                if jid < 10:
                    jfill = "00000" + str(jid)
                elif jid < 100:
                    jfill = "0000" + str(jid)
                elif jid < 1000:
                    jfill = "000" + str(jid)
                elif jid < 10000:
                    jfill = "00" + str(jid)
                elif jid < 100000:
                    jfill = "0" + str(jid)
                else:
                    jfill = str(jid)

                if sampling == "0":

                    jidctrl = 0
                    for l in range(0, lmax + njids):

                        ll = l
                        if ll == 0 or jids[l] == 0 or jidctrl == 1:
                            jids[l] += 1
                        if jids[l] > nruns0[l]:
                            jids[l] = 1
                            jidctrl = 1
                        else:
                            jidctrl = 0

                ipar = []
                if sampling == "0":

                    if float(vhrlmin) == float(vhrlmax):
                        vhrl = float(vhrlmin)
                    else:
                        vhrl = float(vhrlmin) + float(jids[0] - 1) / float(nruns0[0] - 1) * (float(vhrlmax) - float(vhrlmin))
                        ipar.append(1)

                elif int(sampling) > 0:

                    vhrl = round(random.uniform(float(vhrlmin), float(vhrlmax)), 2)

                else:

                    if ltest < lmax + njids:
                        while nvar[ltest] == 0 and ltest < lmax + njids:
                            ltest += 1
                            lnrun = 0

                    if float(vhrlmin) == float(vhrlmax) or not ltest == 0:
                        vhrl = float(vhrl3)
                    else:
                        vhrl = float(vhrlmin) + float(lnrun) / float(nruns00 - 1) * (float(vhrlmax) - float(vhrlmin))

                vhrl = str(vhrl)

                if sampling == "0":

                    if float(vhemmin) == float(vhemmax):
                        vhem = float(vhemmin)
                    else:
                        vhem = float(vhemmin) + float(jids[1] - 1) / float(nruns0[1] - 1) * (float(vhemmax) - float(vhemmin))
                        if phases == "1":
                            ipar.append(2)
                        else:
                            ipar.append(3)

                elif int(sampling) > 0:

                    vhem = round(random.uniform(float(vhemmin), float(vhemmax)), 2)

                else:

                    if float(vhemmin) == float(vhemmax) or not ltest == 1:
                        vhem = float(vhem3)
                    else:
                        vhem = float(vhemmin) + float(lnrun) / float(nruns00 - 1) * (float(vhemmax) - float(vhemmin))

                vhem = str(vhem)

                if ( phases == "3" and hrelease and rhrelease1 ):

                    if sampling == "0":

                        if float(rhrlsmin) == float(rhrlsmax):
                            rhrls = float(rhrlsmin)
                        else:
                            rhrls = float(rhrlsmin) + float(jids[2] - 1) / float(nruns0[2] - 1) * (float(rhrlsmax) - float(rhrlsmin))
                            ipar.append(2)

                    elif int(sampling) > 0:

                        rhrls = random.uniform(float(rhrlsmin), float(rhrlsmax))

                    else:

                        if float(rhrlsmin) == float(rhrlsmax) or not ltest == 2:
                            rhrls = float(rhrls3)
                        else:
                            rhrls = float(rhrlsmin) + float(lnrun) / float(nruns00 - 1) * (float(rhrlsmax) - float(rhrlsmin))

                    rhrls = str(rhrls)

                else:

                    rhrls = "-9999"

                if phases == "3" and hentrmax and rhentrmax1:

                    if sampling == "0":

                        if float(rhemsmin) == float(rhemsmax):
                            rhems = float(rhemsmin)
                        else:
                            rhems = float(rhemsmin) + float(jids[njids - 1] - 1) / float(nruns0[njids - 1] - 1) * (float(rhemsmax) - float(rhemsmin))
                            ipar.append(4)

                    elif int(sampling) > 0:

                        rhems = random.uniform(float(rhemsmin), float(rhemsmax))

                    else:

                        if float(rhemsmin) == float(rhemsmax) or not ltest == njids - 1:
                            rhems = float(rhems3)
                        else:
                            rhems = float(rhemsmin) + float(lnrun) / float(nruns00 - 1) * (float(rhemsmax) - float(rhemsmin))

                    rhems = str(rhems)

                else:

                    rhems = "-9999"

                gt = []

                if sampling == "0":

                    for l in range(0, lmax):
                        if float(flowparam[3 * l]) == float(flowparam[3 * l + 1]):
                            gt.append(float(flowparam[3 * l]))
                        else:
                            gt.append(float(flowparam[3 * l]) + float(jids[l + njids] - 1) / float(nruns0[l + njids] - 1)
                                * (float(flowparam[3 * l + 1]) - float(flowparam[3 * l])))
                            if phases == "1":
                                ipar.append(3 + l)
                            else:
                                ipar.append(5 + l)

                elif int(sampling) > 0:

                    for l in range(0, lmax):
                        gt.append(random.uniform(float(flowparam[2 * l]), float(flowparam[2 * l + 1])))

                else:

                    for l in range(0, lmax):
                        if (float(flowparam[3 * l]) == float(flowparam[3 * l + 1]) or not ltest == l + njids):
                            gt.append(float(flowparam[3 * l + 2]))
                        else:
                            gt.append(float(flowparam[3 * l])+ float(lnrun)/ float(nruns00 - 1)
                                * (float(flowparam[3 * l + 1]) - float(flowparam[3 * l])))

                    lnrun += 1
                    if lnrun == nruns00:
                        lnrun = 0
                        ltest += 1

                writeparam(jid, aflag, eflag, kflag, mflag, tflag, vflag, pf, cores, cellsize, phases, gravity, limiter, layers, ctopo, csurface, centrainment, cstopping, cdispersion, cmelt, aoicoords, 
                    elevation, hrelease, rhrelease1, rhrls, vhrelease, vhrl, hrelease1, hrelease2, hrelease3, trelease, trelstop, stoptime, tslide, zones, impactarea, hdeposit, 
                    vinx, viny, vinx1, viny1, vinx2, viny2, vinx3, viny3, hentrmax, rhentrmax1, rhems, vhentrmax, vhem, hentrmax1, hentrmax2, hentrmax3, 
                    phi1, phi2, phi3, delta1, delta2, delta3, addfri1, addfri2, addfri3, coh1, coh2, coh3, ny1, ny2, ny3, cdeform, zfrag, ambdrag, centr, ctrans12, ctrans13, ctrans23, 
                    ortho_c1, ortho_c2, ortho_c3, ortho_c4, ortho_c5, ortho_c6, hydrograph, hydrocoords, adaptograph, frictiograph, transformograph, gt, sampling, cfl, times, slomo, 
                    cvisc, cflux, thresholds, profile, ctrlpoints, reftime, phexagg, orthophoto, lmax, nmesh, tint, tstop, visualization, ascpath)

                if bingrass0 == "grass":

                    g_region = Module("g.region", n=rnorth, s=rsouth, e=reast, w=rwest, res=cellsize, run_=False)
                    r_avaflow_mult = Module("r.avaflow.mult", jid=str(jid), run_=False)
                    mm = MultiModule([g_region, r_avaflow_mult], sync=False, set_temp_region=True)
                    process_queue.put(mm)

                else:

                    os.mkdir(temppath + "/tmp%s" % jid)
                    strtmp = temppath + "/tmp%s/batch" % jid
                    out = open(strtmp + str(jid), "w")
                    os.system("rm -rf " + locpath + "/map%s" % jid)
                    grass.run_command("g.mapset", flags="c", mapset="map%s" % jid, location=ambvars.LOCATION_NAME)
                    grass.run_command("g.mapsets", mapset=ambvars.MAPSET, operation="add")

                    grass.run_command("g.mapset", mapset=ambvars.MAPSET, location=ambvars.LOCATION_NAME)
                    os.system("mkdir " + locpath + "/map%s/.tmp/rtemp" % jid)
                    os.environ["PATH"] += (os.pathsep + os.path.join(temppath + "/tmp%s") % jid)

                    print(
"""#!/bin/bash
export jid=%s
export SHELL=\"/bin/bash\"
export cellsize=%s
export rnorth=%s
export rsouth=%s
export rwest=%s
export reast=%s
export rtemp=%s
export GRASS_BATCH_JOB=%s/r.avaflow.mult
cp %s/*.txt %s/map%s/.tmp/rtemp/
%s --text %s/map%s
unset GRASS_BATCH_JOB"""
                        % (jid, cellsize, rnorth, rsouth, rwest, reast, temppath, scriptpath2, temppath, locpath, jid, bingrass, locpath, jid), file=out)

                    out.close()
                    fd = os.open(strtmp + str(jid), os.O_RDONLY)
                    os.fchmod(fd, 0o755)
                    os.close(fd)

            start_batch = time.time()

            if bingrass0 == "grass":

                nameList = list(range(1, nruns + 1))
                
                process_queue.wait()
                proc_list = process_queue.get_finished_modules()

                error = 0
                for proc in proc_list:
                    if proc.returncode != 0:
                        grass.error(
                            _("Error running module: %\n    stderr: %s")
                            % (proc.get_bash(), proc.outputs.stderr)
                        )
                        error += 1

                if error > 0:
                    grass.fatal(_("Error running modules."))

                grass.run_command("g.message", message=" ")
                grass.run_command("g.message", message="Batch processing completed.")
                grass.run_command("g.message", message=" ")
                
            else:

                neff = min(ncores, nruns)
                threadList = list(range(1, neff + 1))
                nameList = list(range(1, nruns + 1))
                threads = []
                for tName in threadList:
                    thread = myThread(tName, workQueue)
                    thread.start()
                    threads.append(thread)
                queueLock.acquire()
                for word in nameList:
                    workQueue.put(word)
                queueLock.release()
                while not workQueue.empty():
                    pass
                global exitFlag
                exitFlag = 1
                for t in threads:
                    t.join()
                print()
                print("Batch processing completed.")
                print()

                for jid in range(1, nruns + 1):
                    os.system("rm -rf " + locpath + "/map" + str(jid))

            for j in range(0, 3):
                grass.mapcalc('"%s_iii%s"=0' % (pf, mstring[j]), overwrite=True, quiet=True)
            grass.mapcalc('"%s_dii"=0' % pf, overwrite=True, quiet=True)

            nsuccess = 0
            for jid in nameList:

                ftimesteps = open(pf + "_results/" + pf + "_files/" + pf + "_nout" + str(jid) + ".txt", "r")
                ftimesteps.readline()
                csuccess = (ftimesteps.readline())
                csuccess = int(csuccess.replace("\n", ""))
                ftimesteps.close()

                if jid < 10:
                    jfill = "00000" + str(jid)
                elif jid < 100:
                    jfill = "0000" + str(jid)
                elif jid < 1000:
                    jfill = "000" + str(jid)
                elif jid < 10000:
                    jfill = "00" + str(jid)
                elif jid < 100000:
                    jfill = "0" + str(jid)
                else:
                    jfill = str(jid)

                if csuccess == 1:

                    nsuccess = nsuccess + 1

                if phases == "1":
                    mstringlist = ["_hflow_max", "_hflow_fin", "_vflow_max", "none", "none", "_tflow_max", "_pflow_max", "_basechange_fin", "_treach"]
                elif phases == "3":
                    mstringlist = ["_hflow_max", "_hflow_fin", "_vflow1_max", "_vflow2_max", "_vflow3_max", "_tflow_max", "_pflow_max", "_basechange_fin", "_treach"]

                if csuccess == 1:

                    for mstringi in mstringlist:

                        if not mstringi == "none":

                            grass.run_command("r.in.gdal", input=ascpath + pf + mstringi + str(jid) + ".asc", output=pf + mstringi + str(jid), overwrite=True)

                    minval = grass.raster_info(pf + "_basechange_fin" + str(jid))["min"]
                    maxval = grass.raster_info(pf + "_basechange_fin" + str(jid))["max"]

                    if minval == 0 and maxval == 0:
                        grass.mapcalc('"%s_hflow_dep%s"="%s_hflow_fin%s"' % (pf, str(jid), pf, str(jid)), overwrite=True, quiet=True)
                    else:
                        grass.mapcalc('"%s_hflow_dep%s"=if("%s_basechange_fin%s">0,"%s_basechange_fin%s",0)' % (pf, str(jid), pf, str(jid), pf, str(jid)), overwrite=True, quiet=True)

                    grass.mapcalc('"%s_iis%s"=if("%s_hflow_max%s">%s,1,0)'% (pf, jfill, pf, str(jid), thresholds[0]), overwrite=True, quiet=True)
                    grass.mapcalc('"%s_dis%s"=if("%s_hflow_dep%s">%s,1,0)'% (pf, jfill, pf, str(jid), thresholds[0]), overwrite=True, quiet=True)

                    for j in range(0, 3):

                        grass.mapcalc('"%s_iii%s"=if("%s%s_max%s">%s,"%s_iii%s"+1,"%s_iii%s")' % (pf, mstring[j], pf, mstring[j], str(jid), thresholds[j], pf, mstring[j], pf, mstring[j]),
                            overwrite=True, quiet=True)

                    grass.mapcalc('"%s_dii"=if("%s_hflow_dep%s">%s,"%s_dii"+1,"%s_dii")' % (pf, pf, str(jid), thresholds[0], pf, pf), 
                        overwrite=True, quiet=True)

                    grass.mapcalc('"_hflowii"=if("%s_hflow_max%s">=%s,"%s_hflow_max%s",0)' % (pf, str(jid), thresholds[0], pf, str(jid)), overwrite=True, quiet=True)
                    grass.mapcalc('"_tflowii"=if("%s_tflow_max%s">=%s,"%s_tflow_max%s",0)' % (pf, str(jid), thresholds[1], pf, str(jid)), overwrite=True, quiet=True)

                    grass.run_command("r.out.gdal", input="_hflowii", output=pf + "_results/" + pf + "_aimec/depth/" + jfill + ".asc", format="AAIGrid", overwrite=True)
                    grass.run_command("r.out.gdal", input="_tflowii", output=pf + "_results/" + pf + "_aimec/pressure/" + jfill + ".asc", format="AAIGrid", overwrite=True)

                    corrasc(pf + "_results/" + pf + "_aimec/depth/" + jfill)
                    corrasc(pf + "_results/" + pf + "_aimec/pressure/" + jfill)

                else:
                    grass.mapcalc('"%s_iis%s"=0' % (pf, jfill), overwrite=True, quiet=True)
                    grass.mapcalc('"%s_dis%s"=0' % (pf, jfill), overwrite=True, quiet=True)

            for j in range(0, 3):

                grass.mapcalc('"%s_iii%s"=float("%s_iii%s")/float(%s)' % (pf, mstring[j], pf, mstring[j], str(nsuccess)), overwrite=True, quiet=True)

                grass.run_command("r.out.gdal", input=pf + "_iii" + mstring[j], output=ascpath + pf + "_iii" + mstring[j] + ".asc", format="AAIGrid", 
                    overwrite=True)
                corrasc(ascpath + pf + "_iii" + mstring[j])

            grass.mapcalc('"%s_dii"=float("%s_dii")/float(%s)' % (pf, pf, str(nsuccess)), overwrite=True, quiet=True)

            grass.run_command("r.out.gdal", input=pf + "_dii", output=ascpath + pf + "_dii.asc", format="AAIGrid", overwrite=True)
            corrasc(ascpath + pf + "_dii")

            if impactarea:
                grass.mapcalc('"_cimpactarea"=if("%s">0,1,0)' % impactarea, overwrite=True, quiet=True)
                grass.run_command("r.out.gdal", input="_cimpactarea", output=pf + "_results/" + pf + "_aimec/" + pf + "_impactarea.asc", format="AAIGrid", overwrite=True)
                corrasc(pf + "_results/" + pf + "_aimec/" + pf + "_impactarea")

            if hdeposit:
                grass.mapcalc('"_cdepositarea"=if("%s">0,1,0)' % hdeposit, overwrite=True, quiet=True)
                grass.run_command("r.out.gdal", input="_cdepositarea", output=pf + "_results/" + pf + "_aimec/" + pf + "_depositarea.asc", format="AAIGrid", overwrite=True)
                corrasc(pf + "_results/" + pf + "_aimec/" + pf + "_depositarea")

            stop_batch = time.time()
            comptime_batch = (stop_batch - start_batch)

            timefile = open(filepath + pf + "_time.txt", "w")
            timefile.write(str(comptime_batch))
            timefile.close()

    if vflag:

        if not ( mflag and bingrass0 == "grass"): print()
        print("3. VISUALIZING RESULTS.")
        print()

        ftimesteps = open(pf + "_results/" + pf + "_files/" + pf + "_nout1.txt", "r")
        ntimesteps = ftimesteps.readline()
        ntimesteps = int(ntimesteps.replace("\n", ""))
        csuccess = (ftimesteps.readline())
        csuccess = int(csuccess.replace("\n", ""))
        basechange = (ftimesteps.readline())
        basechange = int(basechange.replace("\n", ""))
       
        ftimesteps.close()

        if mflag:
            success = nsuccess
        else:
            success = csuccess

        if profile and not mflag and success > 0:

            subprocess.call("Rscript %s_results/%s_plots/r.avaflow.profile.R --slave --quiet" % (pf, pf), shell=True)

            for j in [1, 2, 4,]:

                pnames = []
                pnamesc = []
                pwidthc = 400

                for step in range(0, ntimesteps + 1):

                    if step < 10:
                        fill = "000" + str(step)
                    elif step < 100:
                        fill = "00" + str(step)
                    elif step < 1000:
                        fill = "0" + str(step)
                    else:
                        fill = str(step)

                    pimg = Image.open(pf + "_results/" + pf + "_plots/" + pf + "_profiles_timesteps/" + pf + mstring[j] + fill + ".png")
                    pheightc = int(float(pimg.size[1]) * float(pwidthc / float(pimg.size[0])))
                    pimgc = pimg.resize((pwidthc, pheightc), Image.BILINEAR)
                    pimgc.save(temppath + "/" + mstring[j] + "c" + fill + ".png")

                    pnames.append(pf + "_results/" + pf + "_plots/" + pf + "_profiles_timesteps/" + pf + mstring[j] + fill + ".png")
                    pnamesc.append(temppath + "/" + mstring[j] + "c" + fill + ".png")

                images = []

                for nnn in pnames:
                    frame = Image.open(nnn)
                    frame = frame.convert("P", palette=Image.ADAPTIVE, colors=256)
                    images.append(frame)

                images[0].save(pf + "_results/" + pf + "_plots/" + pf + mstring[j] + "_profile.gif", save_all=True, append_images=images[1:], duration=200, loop=0)

                images = []

                for nnn in pnamesc:
                    frame = Image.open(nnn)
                    frame = frame.convert("P", palette=Image.ADAPTIVE, colors=256)
                    images.append(frame)

                images[0].save(pf + "_results/" + pf + "_plots/" + pf + mstring[j] + "_profilec.gif", save_all=True, append_images=images[1:], duration=200, loop=0)

        if (hydrograph or hydrocoords) and not mflag and success > 0:

                subprocess.call("Rscript %s_results/%s_plots/r.avaflow.hydrograph.R --slave --quiet" % (pf, pf), shell=True)

        if success > 0:

            subprocess.call("Rscript %s_results/%s_plots/r.avaflow.map.R --slave --quiet" % (pf, pf), shell=True)

            if not mflag:

                jrange = [0]
                if aflag:
                    jrange.append(1)
                    jrange.append(2)
                if basechange != 0:
                    jrange.append(3)
                if phases == "3" and tflag:
                    jrange.append(5)
                jrange.append(6)

                for j in jrange:

                    mnames = []
                    mnamesc = []
                    mwidthc = 640
            
                    for i in range( 0, ntimesteps + 2 ):

                        if i < 10:
                            fill = "000" + str(i)
                        elif i < 100:
                            fill = "00" + str(i)
                        elif i < 1000:
                            fill = "0" + str(i)
                        else:
                            fill = str(i)

                        if i <= ntimesteps and not j == 6:

                            mimg = Image.open( pf + "_results/" + pf + "_plots/" + pf + "_maps_timesteps/" + pf + mstring[j] + fill + ".png" )
                            mheightc = int(float( mimg.size[1]) * float( mwidthc / float(mimg.size[0] )))
                            mimgc = mimg.resize(( mwidthc, mheightc ), Image.BILINEAR )
                            mimgc.save( temppath + "/" + mstring[j] + "c" + fill + ".png" )
                            mnames.append( pf + "_results/" + pf + "_plots/" + pf + "_maps_timesteps/" + pf + mstring[j] + fill + ".png" )
                            mnamesc.append( temppath + "/" + mstring[j] + "c" + fill + ".png" )

                    if not j == 6:

                        images = []

                        for nnn in mnames:
                            frame = Image.open( nnn )
                            frame = frame.convert( "P", palette=Image.ADAPTIVE, colors=256 )
                            images.append( frame )

                        images[0].save( pf + "_results/" + pf + "_plots/" + pf + mstring[j] + "_map.gif", save_all=True, append_images=images[1:], duration=200, loop=0 )

                        images = []

                        for nnn in mnamesc:
                            frame = Image.open( nnn )
                            frame = frame.convert( "P", palette=Image.ADAPTIVE, colors=256 )
                            images.append( frame )

                        images[0].save( pf + "_results/" + pf + "_plots/" + pf + mstring[j] + "_mapc.gif", save_all=True, append_images=images[1:], duration=200, loop=0 )

            else:

                if impactarea:

                    grass.run_command("r.out.gdal", flags="c", input=impactarea, output=temppath + "/observed.asc", format="AAIGrid", overwrite=True)  
                    grass.run_command("r.out.gdal", flags="c", input=pf + '_iii_hflow', output=temppath + "/xindex.asc", format="AAIGrid", overwrite=True)  

                    with open(temppath + "/xindex.asc", "r") as xfin:
                        with open(temppath + "/index.asc", "w") as xfout:
                            for line in xfin:
                                xfout.write(line.replace("-nan", "0.0"))

                    xfin.close()
                    xfout.close()

                    subprocess.call("Rscript %s_results/%s_plots/r.avaflow.roc.R %s %s %s %s %s --slave --quiet" % (pf, pf, temppath, pf, "1", "1", "1"), shell=True)

                    subprocess.call("Rscript %s_results/%s_plots/r.avaflow.roc.R %s %s %s %s %s --slave --quiet" % (pf, pf, temppath, pf, "2", "2", "1"), shell=True)

                if hdeposit:

                    grass.run_command("r.out.gdal", flags="c", input=hdeposit, output=temppath + "/observed.asc", format="AAIGrid", overwrite=True)  
                    grass.run_command("r.out.gdal", flags="c", input=pf + '_dii', output=temppath + "/xindex.asc", format="AAIGrid", overwrite=True)  

                    with open(temppath + "/xindex.asc", "r") as xfin:
                        with open(temppath + "/index.asc", "w") as xfout:
                            for line in xfin:
                                xfout.write(line.replace("-nan", "0.0"))

                    xfin.close()
                    xfout.close()

                    subprocess.call("Rscript %s_results/%s_plots/r.avaflow.roc.R %s %s %s %s %s --slave --quiet" % (pf, pf, temppath, pf, "3", "1", "2"), shell=True)

                    subprocess.call("Rscript %s_results/%s_plots/r.avaflow.roc.R %s %s %s %s %s --slave --quiet" % (pf, pf, temppath, pf, "4", "2", "2"), shell=True )

                if eflag and mflag and impactarea:

                    if sampling == "0" and len(ipar) == 2:

                        subprocess.call("Rscript %s_results/%s_plots/r.avaflow.multval.R %s %s %s %s %s %s %s --slave --quiet"
                            % (pf, pf, pf, phases, str(nruns), str(ipar[0]), str(ipar[1]), "i", "0"), shell=True)

                if eflag and mflag and hdeposit:

                    if sampling == "0" and len(ipar) == 2:

                        subprocess.call("Rscript %s_results/%s_plots/r.avaflow.multval.R %s %s %s %s %s %s %s --slave --quiet"
                            % (pf, pf, pf, phases, str(nruns), str(ipar[0]), str(ipar[1]), "d", "0"), shell=True)

                if ctrlpoints and reftime:

                    if sampling == "0" and len(ipar) == 2:

                        for ictrlpoint in range(1, len(ctrlpoints) / 2 + 1):

                            subprocess.call("Rscript %s_results/%s_plots/r.avaflow.multval.R %s %s %s %s %s %s %s --slave --quiet" 
                                % (pf, pf, pf, phases, str(nruns), str(ipar[0]), str(ipar[1]), "t", str(ictrlpoint)), shell=True)

    if not vflag:

        ftimesteps = open(pf + "_results/" + pf + "_files/" + pf + "_nout1.txt", "r")
        ftimesteps.readline()
        csuccess = (ftimesteps.readline())
        csuccess = int(csuccess.replace("\n", ""))
        basechange = (ftimesteps.readline())
        basechange = int(basechange.replace("\n", ""))
        ftimesteps.close()

        if mflag:
            success = nsuccess
        else:
            success = csuccess

    os.system("rm -rf " + temppath)
    if eflag and basechange == 0:
        os.system("rm -rf " + ascpath + pf + "_basechange*")

    grass.run_command("g.remove", flags="f", type="rast", pattern="_*", quiet=True)

    if not kflag:
        grass.run_command("g.remove", flags="f", type="rast", pattern=pf + "_*", quiet=True)

    grass.run_command("g.region", flags="d")

    if eflag and not mflag:

        if csuccess == 1:
            print()
            print("Completed successfully in %.1f seconds (net computing time excluding visualization)."% comptime)
            print("Please find the collected results in the directory %s_results." % pf)
            print()

        elif csuccess == 0:
            print()
            print("The simulation was interrupted due to numerical failure.")
            print("Please find the collected results in the directory %s_results." % pf)
            print()

    elif eflag and mflag:

        print()
        print("Completed in %.1f seconds (net computing time excluding visualization)." % comptime_batch)
        print("%s out of %s simulations were successful." % (str(nsuccess), str(nruns)))
        print("Please find the collected results in the directory %s_results." % pf)
        print()

    elif vflag and not eflag:

        print()
        print("Completed successfully.")
        print("Please find the collected results in the directory %s_results." % pf)
        print()

    sys.exit()

if __name__ == "__main__":
    options, flags = grass.parser()
    main()
