#ifndef _write_pdb_atom_record_v2_h
#define _write_pdb_atom_record_v2_h

/* The last two records are compatible with official PBD documentation. 
 * But it appreas that these records can also include charecter id for 
 * segments, i.e. A for protein B for water etc.
 * Character arrays are C style strings. We have to make space for the extra NULL character
 * at the end. So the lenghts of fields are one less than the lengths of the character arrays.*/
typedef struct {
   int      atom_serial_number;
   char     atom_name[5];
   char     atom_alt_loc;
   char     residue_name[4];
   char     chain_id;
   int      residue_sequence_number;
   char     insert_code;
   float    x;
   float    y;
   float    z;
   float    occupancy;
   float    temperature_factor;
   char     element_symbol[3];
   char     charge_on_atom[3];   
}PDB_ATOM;

/* Write one ATOM record to a file.*/
void Write_one_PDB_ATOM(FILE *out, PDB_ATOM *one_PDB_ATOM);



#endif
