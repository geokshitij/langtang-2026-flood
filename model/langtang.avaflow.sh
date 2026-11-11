export AV4_CURVCTRL=3
export AV4_ERO_EXACT=1
export AV4_ERO_BUOY=1
export AV4_ERO_GATE=1
export AV4_ERO_ETA=1
export AV4_ERO_MIXVEL=2

export AV4_ERO_H=50
export AV4_ERO_PREF=0.0275
export AV4_ERO_EMAX=1.0
export AV4_ERO_RHOSB=2200
export AV4_ERO_ALPHASB=0.65
export AV4_ERO_MUSB=0.5
export AV4_ERO_FICE=3
export AV4_ERO_HMIN=2
export AV4_ERO_PSIGN=1
export AV4_ERO_DM=0.5
export AV4_ERO_DB=0.5
export AV4_ERO_DSOLID=2
export AV4_ERO_NPORE=0.3
export AV4_ERO_REBED=1
export AV4_ERO_DSPLIT=1
export AV4_ERO_EBOOK=1
export AV4_ERO_ESPLIT=0
export AV4_ERO_SPEED=1
export AV4_ERO_ALPHASB_RAST=r_alphasb
export AV4_ERO_MUSB_RAST=r_musb

cellsize="30"

g.region -d
g.region res="$cellsize" -a -s

r.in.gdal -o --overwrite --quiet input=inputs/maps/hotstart_bed_dem.tif output=r_elev0
r.in.gdal -o --overwrite --quiet input=inputs/maps/release_rock_utm45n.tif output=r_hrock0
r.in.gdal -o --overwrite --quiet input=inputs/maps/release_ice_utm45n.tif output=r_hice0
r.in.gdal -o --overwrite --quiet input=inputs/maps/hotstart_water_depth_n025.tif output=r_hwater0
r.in.gdal -o --overwrite --quiet input=inputs/maps/hotstart_water_vx_n025.tif output=r_vinx30
r.in.gdal -o --overwrite --quiet input=inputs/maps/hotstart_water_vy_n025.tif output=r_viny30
r.in.gdal -o --overwrite --quiet input=inputs/maps/entrrock_n40.tif output=r_entrrock0
r.in.gdal -o --overwrite --quiet input=inputs/maps/glacier_ice_utm45n.tif output=r_entrice0
r.in.gdal -o --overwrite --quiet input=inputs/maps/entrwater_n40.tif output=r_entrwater0
r.in.gdal -o --overwrite --quiet input=inputs/maps/fragmentation_zone_utm45n.tif output=r_zfrag0
r.in.gdal -o --overwrite --quiet input=inputs/maps/deltaz_T3.tif output=r_delta0
r.in.gdal -o --overwrite --quiet input=inputs/maps/addfriz_T3.tif output=r_addfri0
r.in.gdal -o --overwrite --quiet input=inputs/maps/alphasb_T3.tif output=r_alphasb
r.in.gdal -o --overwrite --quiet input=inputs/maps/musb_T3.tif output=r_musb

r.mapcalc --overwrite "r_elev = if(isnull(r_elev0), 0, r_elev0)"
r.mapcalc --overwrite "r_hrock = if(isnull(r_hrock0), 0, r_hrock0)"
r.mapcalc --overwrite "r_hice = if(r_hrock > 0, if(isnull(r_hice0), 0, r_hice0), 0)"
r.mapcalc --overwrite "r_hwater = if(isnull(r_hwater0), 0, r_hwater0)"
r.mapcalc --overwrite "r_vinx3 = if(isnull(r_vinx30), 0, r_vinx30)"
r.mapcalc --overwrite "r_viny3 = if(isnull(r_viny30), 0, r_viny30)"
r.mapcalc --overwrite "r_entrrock = if(isnull(r_entrrock0), 0, r_entrrock0)"
r.mapcalc --overwrite "r_entrice = if(r_hrock + r_hice > 0, 0, if(isnull(r_entrice0), 0, r_entrice0))"
r.mapcalc --overwrite "r_entrwater = if(isnull(r_entrwater0), 0, r_entrwater0)"
r.mapcalc --overwrite "r_zfrag = if(isnull(r_zfrag0), 0, r_zfrag0)"
r.mapcalc --overwrite "r_delta = if(isnull(r_delta0), 7, r_delta0)"
r.mapcalc --overwrite "r_addfri = if(isnull(r_addfri0), -2000, r_addfri0)"

hydrographs="inputs/hydrographs/AUG2_T00.txt,inputs/hydrographs/AUG2_T01.txt,inputs/hydrographs/AUG2_T02.txt,inputs/hydrographs/AUG2_T03.txt,inputs/hydrographs/AUG2_T04.txt,inputs/hydrographs/AUG2_T05.txt,inputs/hydrographs/AUG2_T06.txt,inputs/hydrographs/AUG2_T07.txt,inputs/hydrographs/AUG2_T08.txt,inputs/hydrographs/AUG2_T09.txt,inputs/hydrographs/AUG2_T10.txt,inputs/hydrographs/AUG2_T11.txt,inputs/hydrographs/AUG2_T12.txt,inputs/hydrographs/AUG2_T13.txt,inputs/hydrographs/AUG2_T14.txt,inputs/hydrographs/AUG2_T15.txt,inputs/hydrographs/AUG2_T16.txt"
hydrocoords="353745.0,3131055.0,1,-1,351945.0,3134835.0,1,-1,347355.0,3135705.0,1,-1,346485.0,3134535.0,1,-1,340665.0,3129525.0,1,-1,336855.0,3118485.0,1,-1,337575.0,3116325.0,1,-1,333135.0,3112125.0,1,-1,323655.0,3106815.0,1,-1,320565.0,3096645.0,1,-1,321735.0,3095865.0,1,-1,316605.0,3089775.0,1,-1,314025.0,3082995.0,1,-1,309705.0,3083205.0,1,-1,304875.0,3080385.0,1,-1,304785.0,3077835.0,1,-1,303225.0,3076485.0,1,-1"

r.avaflow.40G -ea prefix=FINAL cellsize="$cellsize" phases=3 elevation=r_elev hrelease1=r_hrock hrelease2=r_hice hrelease3=r_hwater vinx3=r_vinx3 viny3=r_viny3 hentrmax1=r_entrrock hentrmax2=r_entrice hentrmax3=r_entrwater hydrograph="$hydrographs" hydrocoords="$hydrocoords" zfrag=r_zfrag delta1=r_delta delta2=r_delta addfri1=r_addfri addfri2=r_addfri centrainment=-1 density=2700,1000,1080 friction=40,40,0,10,10,0,-2000,-2000,0.025 viscosity=-7,-7,-3 deformation=0.5,0.5,1.0 fragmentation=3,0 thresholds=1,100000,100000,1,0.001 time="${TIME:-30,24000}" cmelt=1 melting=-15,+5,+5,0.0025,0.2

g.region -d
