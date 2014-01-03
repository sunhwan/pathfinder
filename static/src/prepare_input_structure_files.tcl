# Last modified on December 31, 2013.
# This is an updated version that allows more flexibility.
#
# Prepares two PDBs with CAs from user input.
# This script is intended for VMD.
#
#
# User inputs are stored in two files, one for each end state, called 
# 'INPUT_STRUCTURE_INFO_ENDPOINT_1' and 'INPUT_STRUCTURE_INFO_ENDPOINT_2'.
# Format:
# -----------------------------------------------------------------
# (Input PDB file for end point)
# (Chain id) (Residue numbers to include for this chain, VMD syntax) [One chain per line]
# -----------------------------------------------------------------
#
# Tasks performed by this script
# 1. Read the 'INPUT_STRUCTURE_INFO_ENDPOINT_1' and 'INPUT_STRUCTURE_INFO_ENDPOINT_2' files
# 2. Extract C-alphas for the given selections and align them with C-alphas of structure 1 as reference
# 3. Write two PDB files with aligned coordinates called '1_CA.pdb' and '2_CA.pdb'.
# 4. Write two myxyz files called 'INPUT_STRUCTURE_1' and 'INPUT_STRUCTURE_2'. 
#    Also write the 'REFERENCE_FOR_ALIGNMENT' file which is same as 'INPUT_STRUCTURE_1'.
# 5. Write the 'PDB_INFO' file from '1_CA.pdb'
# 6. Write a log file with some output.
###################################################################################################################


########## Tcl procedures ###########

# Procedure for writing a myxyz file from a VMD atom selection
proc Write_myxyz_from_vmd_atom_selection {sel filename_myxyz} {
    set coords [$sel get {x y z}]
    set out [open $filename_myxyz w]
    foreach atom_coord $coords {
        puts $out $atom_coord
    }
    close $out
}


# Write the PDB_INFO file
# Format:
# ---------------------------------
#  Number of atoms (%d)
#  Each line will have the follwoing fields separated by white space
#  (Atom_serial_number) (Atom_name) (Three_letter_residue_name) (Chain_identifier) (Residue_serial_number) (Element_symbol)
# ---------------------------------
proc Write_PDB_INFO {file_name_PDB} {
    mol load pdb $file_name_PDB
    set sel [atomselect top "all"]

    set NUM_ATOMS [$sel num]

    set SERIAL_NUMBERS [$sel get serial]
    set ATOM_NAMES [$sel get name]
    set RESIDUE_NAMES [$sel get resname]
    set CHAIN_IDS [$sel get chain]
    set RESIDUE_NUMBERS [$sel get resid]
    set ELEMENT_SYMBOLS [$sel get element]

    set out [open "PDB_INFO" w]
    puts $out $NUM_ATOMS
    for {set i 0} {$i < $NUM_ATOMS} {incr i} {
        puts $out "[lindex $SERIAL_NUMBERS $i] [lindex $ATOM_NAMES $i] [lindex $RESIDUE_NAMES $i] [lindex $CHAIN_IDS $i] [lindex $RESIDUE_NUMBERS $i] [lindex $ELEMENT_SYMBOLS $i]"
    }
    close $out
}

# Read 'INPUT_STRUCTURE_INFO_ENDPOINT_1' or 'INPUT_STRUCTURE_INFO_ENDPOINT_2'
# Populates through upvar the following entires: Name of the PDB file
#                                                Chain IDs as a list
#                                                Residue ranges of all chain ids as a list
proc Read_strucutre_info_file_for_one_end_point { structure_info_filename input_pdb_filename chain_id_list residue_selection_list  } {
    upvar $input_pdb_filename input_pdb_filename_local
    upvar $chain_id_list chain_id_list_local
    upvar $residue_selection_list residue_selection_list_local

    set in [open $structure_info_filename r]
    set file_data [read $in]
    close $in

    # Extract all lines
    set all_lines [split $file_data "\n"]
    set num_elems [llength $all_lines]
    set num_lines [expr {$num_elems - 1}]
    # Loop over all the lines
    for {set i 0} {$i < $num_lines} {incr i} {
       set temp_list ""
       set one_line [lindex $all_lines $i]
       # The RE matches any non-empty sequence of non-whitespace characters
       set one_line_splitted [regexp -all -inline {\S+} $one_line]
       if {$i == 0} {
          set input_pdb_filename_local [lindex $one_line_splitted 0]
       } else {
          set num_elements_one_line [llength $one_line_splitted]
          lappend chain_id_list_local [lindex $one_line_splitted 0]
          for {set j 1} {$j < $num_elements_one_line} {incr j} {
              lappend temp_list [lindex $one_line_splitted $j]
          }
          lappend residue_selection_list_local $temp_list
       }
    }
    
}


# Construct the VMD selection string from chain ids and residue ranges
proc Construct_atom_selection_string { num_chains chain_id_list residue_selection_list final_selection  } {
    upvar $final_selection final_selection_local
    set final_selection_local ""
    for {set i 0} {$i < $num_chains} {incr i} {
        set final_selection_local [concat $final_selection_local "(chain " [lindex $chain_id_list $i] " and name CA and resid " [lindex $residue_selection_list $i] ") "]
        if {$i != [expr {$num_chains - 1}]} {
           set final_selection_local [concat $final_selection_local " or "]
        }
    }
}

########## End of Tcl procedures ##########

###################################################################################################################

# Read 'INPUT_STRUCTURE_INFO_ENDPOINT_1'.
# Extract PDB file name, chain ids and atom selections.
set structure_info_filename_1 INPUT_STRUCTURE_INFO_ENDPOINT_1
Read_strucutre_info_file_for_one_end_point $structure_info_filename_1 input_pdb_filename_1 chain_id_list_1 residue_selection_list_1

# Read 'INPUT_STRUCTURE_INFO_ENDPOINT_2'.
# Extract PDB file name, chain ids and atom selections.
set structure_info_filename_2 INPUT_STRUCTURE_INFO_ENDPOINT_2
Read_strucutre_info_file_for_one_end_point $structure_info_filename_2 input_pdb_filename_2 chain_id_list_2 residue_selection_list_2

# Check number of chanins in both selection.
set num_chains_1 [llength $chain_id_list_1]
set num_chains_2 [llength $chain_id_list_2]
if { $num_chains_1 != $num_chains_2 } {
   puts "ERROR: Same number of chains should be selected from both PDB files"
   exit
} else {
   set num_chains $num_chains_1
}


# Construct the VMD atom seelction strings for both PDB files.
Construct_atom_selection_string $num_chains $chain_id_list_1 $residue_selection_list_1 final_selection_string_1
Construct_atom_selection_string $num_chains $chain_id_list_2 $residue_selection_list_2 final_selection_string_2


# Read two input PDB files.
mol load pdb $input_pdb_filename_1
set sel_1 [atomselect top "$final_selection_string_1"]
mol load pdb $input_pdb_filename_2
set sel_2 [atomselect top "$final_selection_string_2"]


# Number of atoms in two selections
set num_atoms_1 [$sel_1 num]
set num_atoms_2 [$sel_2 num]
if {$num_atoms_1 != $num_atoms_2} {
    puts "ERROR: Different number of atoms were extracted from two structures based on input atom selections. Make sure both selections result in same number of atoms."
    exit
} else {
    set num_atoms $num_atoms_1
}


# Calculate the initial RMSD without alignment.
set initial_rmsd [measure rmsd $sel_1 $sel_2]


# Align
set transformation_matrix [measure fit $sel_2 $sel_1]
$sel_2 move $transformation_matrix


# Calculate the final RMSD after alignment
set rmsd_after_alignment [measure rmsd $sel_1 $sel_2]


# Write two PDB files with C-alphas
$sel_1 writepdb 1_CA.pdb
$sel_2 writepdb 2_CA.pdb


# Write myxyz files: 'INPUT_STRUCTURE_1', 'INPUT_STRUCTURE_1' and 'REFERENCE_FOR_ALIGNMENT' (which is same as 'INPUT_STRUCTURE_1')
set filename_myxyz INPUT_STRUCTURE_1
Write_myxyz_from_vmd_atom_selection $sel_1 $filename_myxyz
set filename_myxyz INPUT_STRUCTURE_2
Write_myxyz_from_vmd_atom_selection $sel_2 $filename_myxyz
set filename_myxyz REFERENCE_FOR_ALIGNMENT
Write_myxyz_from_vmd_atom_selection $sel_1 $filename_myxyz


# Write the 'PDB_INFO' file
set file_name_PDB 1_CA.pdb
Write_PDB_INFO $file_name_PDB


# Write a log file called 'structure_prep.log'
set out [open "structure_prep.log" w]
puts $out "Version 2"
puts $out "Input filename for end point 1: $input_pdb_filename_1"
puts $out "Input filename for end point 2: $input_pdb_filename_2"
puts $out "Number of chains selected: $num_chains"
puts $out "String used for VMD atom selection for endpoint 1: $final_selection_string_1"
puts $out "String used for VMD atom selection for endpoint 2: $final_selection_string_2"
puts $out "Number of atoms: $num_atoms"
puts $out "Intial RMSD: $initial_rmsd"
puts $out "RMSD after alignment: $rmsd_after_alignment"
close $out

exit











 
