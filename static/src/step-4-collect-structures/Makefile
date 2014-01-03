#Makefile 'collect_calc_energy_v2.c'
#Last updated May 20, 2013

#GSLINC = /usr/local/include
#GSLPATH = /usr/local/lib

CC = gcc
CFLAGS = -Wall -Wno-unused-result -O3
LDFLAGS = -lm



PROG = 4_collec_ener_v2
HDRS = calc_energy_pathway_v2.h calc_rmsd.h collect_structures.h force_pot_simple_enm.h include_all.h md_ld_common_d_dimension_namd.h read_configuration.h wrtie_configuration.h
SRCS = collect_calc_energy_v2.c calc_energy_pathway_v2.c calc_rmsd.c collect_structures.c force_pot_simple_enm.c read_configuration.c wrtie_configuration.c


OBJS = $(SRCS:.c=.o) 


#Build the executable
$(PROG) : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(PROG)

calc_energy_pathway_v2.o : calc_energy_pathway_v2.h

calc_rmsd.o : calc_rmsd.h

collect_structures.o : collect_structures.h

force_pot_simple_enm.o : force_pot_simple_enm.h

read_configuration.o : read_configuration.h

wrtie_configuration.o : wrtie_configuration.h

collect_calc_energy_v2.o : include_all.h

