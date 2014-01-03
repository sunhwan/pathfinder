# Make file for 'minimize_on_cusp_v4.c'

#Last updated May 20, 2013

#GSLINC = /usr/local/include
#GSLPATH = /usr/local/lib

CC = gcc
CFLAGS = -Wall -Wno-unused-result -O3
LDFLAGS = -lm



PROG = 2_find_min_on_cusp_v4
HDRS = align_two_structures_Lei.h calc_rmsd.h check_step_sizes.h force_pot_simple_enm.h include_all.h md_ld_common_d_dimension_namd.h membership_and_locate_point_on_cusp_v1.h read_configuration.h wrtie_configuration.h 
SRCS = minimize_on_cusp_v4.c align_two_structures_Lei.c calc_rmsd.c check_step_sizes.c force_pot_simple_enm.c membership_and_locate_point_on_cusp_v1.c read_configuration.c wrtie_configuration.c


OBJS = $(SRCS:.c=.o) 


#Build the executable
$(PROG) : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(PROG)

align_two_structures_Lei.o : align_two_structures_Lei.h

calc_rmsd.o : calc_rmsd.h

check_step_sizes.o : check_step_sizes.h

force_pot_simple_enm.o : force_pot_simple_enm.h

membership_and_locate_point_on_cusp_v1.o : membership_and_locate_point_on_cusp_v1.h

read_configuration.o : read_configuration.h

wrtie_configuration.o : wrtie_configuration.h

minimize_on_cusp_v4.c : include_all.h



