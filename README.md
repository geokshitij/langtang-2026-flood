# The 2026 Langtang rock-ice avalanche flood: code and inputs

> **Note.** The paper that this code belongs to is under review. Please do not use, cite or share this code or its
> results before the paper is published.

Code and inputs of the simulation in "Catastrophic cascade of the 2026 Langtang rock-ice avalanche flood, Nepal"
by Dahal et al. The simulation follows the cascade from the collapse on Langtang Lirung to Galchhi, 101 km
downstream.

| Path | Content |
|---|---|
| `model/openmp/`, `model/serial/` | r.avaflow 4.0 with the erosion code and our changes (parallel and serial builds, identical results) |
| `model/langtang.avaflow.sh` | the run, with every setting |
| `model/inputs/` | 14 maps and 17 inflow hydrographs, with checksums |
| `analysis/places.py`, `analysis/data/path.csv` | the hazard numbers of the paper from the run output |

## 1. System requirements

* Linux on x86-64. Tested on Rocky Linux 9.8 with GRASS GIS 8.4.3 in the Docker image
  `osgeo/grass-gis:releasebranch_8_4-ubuntu` (Ubuntu 24.04, gcc 13.3, R 4.3.3), run through Apptainer.
* GRASS GIS 8.4, gcc 12 or newer, make and R (r.avaflow calls R).
* Python 3 with NumPy and SciPy for `places.py`. Tested with Python 3.11.13, NumPy 2.4.6 and SciPy 1.17.1.
* No special hardware. The run needs 18 GB of memory. The demo writes 7 GB, the full run about 800 GB.

## 2. Installation guide

1. Install GRASS GIS 8.4 (https://grass.osgeo.org/download/) or use the Docker image above.
2. In GRASS, install the solver from one of the two folders:

   ```bash
   g.extension extension=r.avaflow url=/path/to/langtang-2026-flood/model/openmp
   ```

The solver compiles in about 15 s. With GRASS GIS installed, the installation takes under 1 min.

## 3. Demo

The first 30 s of the simulation, outputs every 5 s, from `model/`:

```bash
grass -c inputs/maps/hotstart_bed_dem.tif -e grassdata/langtang
TIME=5,30 OMP_NUM_THREADS=16 grass grassdata/langtang/PERMANENT --exec bash langtang.avaflow.sh
```

Expected output: the folder `FINAL_results`, with the flow maps in `FINAL_ascii/` and the tables in `FINAL_files/`.
The last row of `FINAL_files/FINAL_volumes.txt` (t = 30.1 s) holds the flowing rock, ice and water:
106,735,480, 10,173,181 and 11,774,759 m3.

Run time: 8 to 13 min on 16 cores with `openmp/`, 30 min on one core with `serial/`.

## 4. Instructions for use

To run on other data, replace the maps in `model/inputs/` and the settings in `langtang.avaflow.sh`. The r.avaflow
manual (https://www.landslidemodels.org/r.avaflow/) describes every parameter.

To compute the numbers for another run or other places:

```bash
python analysis/places.py RESULTS_DIR PREFIX --tint 30 --places 20.49:34.41
```

`--places` takes distances along the river path in km (`analysis/data/path.csv`).

## Reproduction

1. Run the full simulation from `model/`, in the GRASS location of the demo, 400 min of flow, outputs every 30 s
   (17.6 h on 32 cores):

   ```bash
   OMP_NUM_THREADS=32 grass grassdata/langtang/PERMANENT --exec bash langtang.avaflow.sh
   ```

2. Compute the numbers of the paper:

   ```bash
   python analysis/places.py model/FINAL_results FINAL --out places
   ```

   This writes `places.csv` and `places_series.npz` and prints the table.

## Changes to r.avaflow 4.0

Fixes, always on:
* `pow(x, 1/2)` and `pow(x, 1/3)` in the erosion code (ten places) are square and cube roots.
* The hydrograph start times are no longer written past the end of their array.
* A line of the parameter file may be longer than 1010 characters.

Corrections, each off by default and set in `langtang.avaflow.sh`:

| Variable | Run value | Change |
|---|---|---|
| `AV4_CURVCTRL` | 3 | bed curvature with the physical sign, centripetal load along the flow (Fischer et al. 2012) |
| `AV4_ERO_EXACT` | 1 | the bed loses exactly what the flow gains |
| `AV4_ERO_BUOY` | 1 | one buoyancy factor in the erosion stresses (Pudasaini 2025b, eqs. 59, 61) |
| `AV4_ERO_GATE` | 1 | the erosion depth limit applies to the solid depth |
| `AV4_ERO_ETA` | 1 | dynamic viscosity in the fluid drift (eq. 35) |
| `AV4_ERO_MIXVEL` | 2 | eroded material joins at the volume-weighted velocity (eqs. 46, 47) |

## Licence and credits

GNU General Public License, version 3 (`LICENSE`). r.avaflow 4.0 is by Martin Mergili and Shiva P. Pudasaini with the
contributors listed in the source; the erosion code is Pudasaini's. Please cite the paper, r.avaflow (Mergili et al.
2025) and the models it implements: Pudasaini and Mergili (2019), Pudasaini et al. (2024), Pudasaini (2025a) and
Pudasaini (2025b).
