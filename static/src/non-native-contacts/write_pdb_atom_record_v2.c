/*Codes for writing PDB ATOM record and PDB file given the coordinates and other relevant information.*/

/*Last modified on July 06, 2012*/


/* FORTRAN format for writing PDB records
 * ATOM
 * HETATM Format(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)
 *
 * HELIX  Format(A6,1X,I3,1X,A3,2(1X,A3,1X,A1,1X,I4,A1),I2,A30,1X,I5)
 *
 * SHEET Format(A6,1X,I3,1X,A3,I2,2(1X,A3,1X,A1,I4,A1),I2,2(1X,A4,A3,1X,A1,I4,A1))
 *
 * SSBOND Format(A6,1X,I3,1X,A3,1X,A1,1X,I4,A1,3X,A3,1X,A1,1X,I4,A1,23X,2(2I3,1X),F5.2)
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "write_pdb_atom_record_v2.h"


/* THIS HAS ALL THE FIELDS. Write one ATOM record to a file.
 * FORTRAM Format(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)
 * The code follows this FORTRAN format which is consistent with the official format.
 * */
void Write_one_PDB_ATOM(FILE *out, PDB_ATOM *one_PDB_ATOM)
{
  fprintf(out, "ATOM  %5d %-4s%c%-3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %-2s%-2s\n", 
               one_PDB_ATOM->atom_serial_number, 
               one_PDB_ATOM->atom_name, 
               one_PDB_ATOM->atom_alt_loc, 
               one_PDB_ATOM->residue_name, 
               one_PDB_ATOM->chain_id, 
               one_PDB_ATOM->residue_sequence_number, 
               one_PDB_ATOM->insert_code, 
               one_PDB_ATOM->x, 
               one_PDB_ATOM->y, 
               one_PDB_ATOM->z, 
               one_PDB_ATOM->occupancy, 
               one_PDB_ATOM->temperature_factor, 
               one_PDB_ATOM->element_symbol, 
               one_PDB_ATOM->charge_on_atom);
}




