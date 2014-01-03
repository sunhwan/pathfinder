# Makefile for 'find_close_pair_pathway_v2.c'.
#
# Last updated on August 23, 2013

#GSLINC = /usr/local/include
#GSLPATH = /usr/local/lib

CC = gcc
CFLAGS = -Wall -Wno-unused-result -O3
LDFLAGS = -lm
#CFLAGS = -Wall -O2 -O3 -fopenmp -I$(GSLINC)
#LDFLAGS = -L$(GSLPATH) -lgsl -lgomp -lm


PROG = find_pairs_path_v2
HDRS = include_all.h md_ld_common_d_dimension_namd.h read_configuration.h write_pdb_atom_record_v2.h write_PDB_files.h
SRCS = find_close_pair_pathway_v2.c read_configuration.c write_pdb_atom_record_v2.c write_PDB_files.c


OBJS = $(SRCS:.c=.o) 


#Build the executable
$(PROG) : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(PROG)

read_configuration.o : read_configuration.h

write_pdb_atom_record_v2.o : write_pdb_atom_record_v2.h

write_PDB_files.o : write_PDB_files.h

find_close_pair_pathway_v2.o : include_all.h
