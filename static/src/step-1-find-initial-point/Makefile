# Make file for 'find_structure_on_cusp_v2.c'

#Last updated May 20, 2013

#GSLINC = /usr/local/include
#GSLPATH = /usr/local/lib

CC = gcc
CFLAGS = -Wall -Wno-unused-result -O3
LDFLAGS = -lm



PROG = 1_locate_struct_on_cusp_v2
HDRS = force_pot_simple_enm.h include_all.h md_ld_common_d_dimension_namd.h membership_and_locate_point_on_cusp_v1.h read_configuration.h wrtie_configuration.h 
SRCS = find_structure_on_cusp_v2.c force_pot_simple_enm.c membership_and_locate_point_on_cusp_v1.c read_configuration.c wrtie_configuration.c


OBJS = $(SRCS:.c=.o) 


#Build the executable
$(PROG) : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(PROG)

force_pot_simple_enm.o : force_pot_simple_enm.h

membership_and_locate_point_on_cusp_v1.o : membership_and_locate_point_on_cusp_v1.h

read_configuration.o : read_configuration.h

wrtie_configuration.o : wrtie_configuration.h

find_structure_on_cusp_v2.o : include_all.h



