# Makefile for 'make_pathway_PDB_v2.c'.
#
# Last updated on June 12, 2013

#GSLINC = /usr/local/include
#GSLPATH = /usr/local/lib

CC = gcc
CFLAGS = -Wall -Wno-unused-result -O3
LDFLAGS = -lm
#CFLAGS = -Wall -O2 -O3 -fopenmp -I$(GSLINC)
#LDFLAGS = -L$(GSLPATH) -lgsl -lgomp -lm


PROG = 5_makepathway_v2
HDRS = align_two_structures_using_subset_Lei.h include_all.h md_ld_common_d_dimension_namd.h read_configuration.h write_images_to_PDB_traj_v2.h write_pdb_atom_record_v2.h write_PDB_files.h wrtie_configuration.h
SRCS = make_pathway_PDB_v2.c align_two_structures_using_subset_Lei.c read_configuration.c write_images_to_PDB_traj_v2.c write_pdb_atom_record_v2.c write_PDB_files.c wrtie_configuration.c


OBJS = $(SRCS:.c=.o) 


#Build the executable
$(PROG) : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(PROG)

align_two_structures_using_subset_Lei.o : align_two_structures_using_subset_Lei.h

read_configuration.o : read_configuration.h

write_images_to_PDB_traj_v2.o : write_images_to_PDB_traj_v2.h

write_pdb_atom_record_v2.o : write_pdb_atom_record_v2.h

write_PDB_files.o : write_PDB_files.h

wrtie_configuration.o : wrtie_configuration.h

make_pathway_PDB_v2.o : include_all.h
