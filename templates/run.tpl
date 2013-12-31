#!/bin/csh

set vmd = {{ vmd }}

$vmd -dispdev none -e prepare_input_structure_files.tcl > prepare_input_structure_files.out

set np = `sed -n 's/Number of atoms: \([0-9]*\)/\1/p' structure_prep.log`
set fc1 = {{ fc1 }}
set cutoff1 = {{ cutoff1 }}
set offset1 = {{ offset1 }}
set fc2 = {{ fc2 }}
set cutoff2 = {{ cutoff2 }}
set offset2 = {{ offset2 }}
set rmsd = {{ target }}
set stepsize_cusp1 = {{ stepsize_cusp1 }}
set stepsize_cusp2 = {{ stepsize_cusp2 }}
set stepsize_slid1 = {{ stepsize_slid1 }}
set stepsize_slid2 = {{ stepsize_slid2 }}

perl -pi -e 's/\$FORCE_CONSTANT_1 = ([\.0-9]+)?;/\$FORCE_CONSTANT_1 = '$fc1';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$CUT_OFF_1 = ([\.0-9]+)?;/\$CUT_OFF_1 = '$cutoff1';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$ENERGY_OFFSET_1 = ([\.0-9]+)?;/\$ENERGY_OFFSET_1 = '$offset1';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$FORCE_CONSTANT_2 = ([\.0-9]+)?;/\$FORCE_CONSTANT_2 = '$fc2';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$CUT_OFF_2 = ([\.0-9]+)?;/\$CUT_OFF_2 = '$cutoff2';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$ENERGY_OFFSET_2 = ([\.0-9]+)?;/\$ENERGY_OFFSET_2 = '$offset2';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$NUM_PARTICLES = ([0-9]+)?;/\$NUM_PARTICLES = '$np';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$RMSD_PATHWAY = ([\.0-9]+)?;/\$RMSD_PATHWAY = '$rmsd';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$STEP_SIZE_MINIMIZE_ON_CUSP_1 = ([\.0-9]+)?;/\$STEP_SIZE_MINIMIZE_ON_CUSP_1 = '$stepsize_cusp1';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$STEP_SIZE_MINIMIZE_ON_CUSP_2 = ([\.0-9]+)?;/\$STEP_SIZE_MINIMIZE_ON_CUSP_2 = '$stepsize_cusp2';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$STEP_SIZE_SLIDE_DOWN_1 = ([\.0-9]+)?;/\$STEP_SIZE_SLIDE_DOWN_1 = '$stepsize_slid1';/g' create_input_files_ANMPathway.pl
perl -pi -e 's/\$STEP_SIZE_SLIDE_DOWN_2 = ([\.0-9]+)?;/\$STEP_SIZE_SLIDE_DOWN_2 = '$stepsize_slid2';/g' create_input_files_ANMPathway.pl
./create_input_files_ANMPathway.pl

./1_locate_struct_on_cusp_v2

./2_find_min_on_cusp_v4

cp INPUT_SLIDE_ONE_SURFACE_1 INPUT_SLIDE_ONE_SURFACE
./3_desc_one_surface_v3

rm -f INPUT_SLIDE_ONE_SURFACE
cp INPUT_SLIDE_ONE_SURFACE_2 INPUT_SLIDE_ONE_SURFACE
./3_desc_one_surface_v3

./4_collec_ener_v2

set ni = `sed -n 's/Total number of structures in the pathway: \([0-9]*\)/\1/p' number_of_structures_in_pathway`
echo $np $ni
./5_makepathway_v2 -n $np -i $ni -a 0

echo "NUM_PARTICLES $np\
NUM_IMAGES $ni\
CUT_OFF 5.0\
SEPERATION_END_STATES 10.0" > INPUT_FIND_CLOSE_PAIRS
./find_pairs_path_v2

$vmd -dispdev none -e movie.tcl > movie.out
convert -loop 4 frame.*.rgb movie.gif

tar -czf COORDS_IMAGES.tar.gz COORDS_IMAGE_*
tar -czf COORDS_SURFACE.tar.gz OUT_COORDS_SURFACE_*
tar -czf contact_pairs.tar.gz distance_*
rm -f COORDS_IMAGE_*
rm -f OUT_COORDS_SURFACE_*
rm -f distance_*
rm -f *.rgb

