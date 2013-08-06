#!/bin/csh

set vmd = {{ vmd }}

$vmd -dispdev none -e prepare_input_strucutre_files.tcl > prepare_input_structure_files.out
set np = `sed -n 's/Number of atoms: \([0-9]*\)/\1/p' structure_prep.log`

perl -pi -e 's/\$NUM_PARTICLES = [0-9]+/\$NUM_PARTICLES = '$np'/g' create_input_files_ANMPathway.pl
./create_input_files_ANMPathway.pl

cat '1' > job.progress
./1_locate_struct_on_cusp_v2

cat '2' > job.progress
./2_find_min_on_cusp_v4

cat '3.1' > job.progress
cp INPUT_SLIDE_ONE_SURFACE_1 INPUT_SLIDE_ONE_SURFACE
./3_desc_one_surface_v3

cat '3.2' > job.progress
cp INPUT_SLIDE_ONE_SURFACE_2 INPUT_SLIDE_ONE_SURFACE
./3_desc_one_surface_v3

cat '4' > job.progress
./4_collec_ener_v2

cat '5' > job.progress
set ni = `sed -n 's/Total number of structures in the pathway: \([0-9]*\)/\1/p' number_of_structures_in_pathway`
echo $np $ni
./5_makepathway_v2 -n $np -i $ni -a 0

tar -cvf COORDS_IMAGES.tar.gz COORDS_IMAGE_*
tar -cvf COORDS_SURFACE.tar.gz OUT_COORDS_SURFACE_*
rm -f COORDS_IMAGE_*
rm -f OUT_COORDS_SURFACE_*