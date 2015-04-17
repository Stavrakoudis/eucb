##############################################################
### Makefile for the program #################################
include ../Makefile
INC = -I../include
OBJS = rad.o command.o command_sasa.o command_writedcd.o command_clusterrmsd.o command_distance.o command_bturn.o command_complexsalt.o command_rmsf.o command_angle.o command_noe.o command_stack.o command_salt.o command_sidedist.o command_anal.o command_backbone.o command_donor.o command_watbridge.o command_unary.o command_aronh.o command_aropos.o command_aroch.o psf.o command_rmsd.o command_center.o command_pdbwrite.o command_pdbangle.o command_pdo.o command_contact.o command_hydrosolve.o command_closewater.o dcd.o team.o rmsd.o options.o donors.o util.o global.o getoptions.o main.o helpfiles.o pdb.o tex.o counter.o command_iwcn.o command_ahelix.o command_watpdb.o eigen.o command_entropy.o parfile.o command_energy.o command_disu.o command_hbonds1.o

alldist: filereader $(OBJS)  $(PROGRAMNAME) 

filereader: filereader.cc
	$(CXX) filereader.cc -o filereader

%.o: %.cc
	$(CXX) $(CFLAGS) $(INC) -c $< -o $@

$(PROGRAMNAME):  $(OBJS)
	./filereader
	$(CXX) $(INC) $(OBJS)  -o $(PROGRAMNAME)
	mv $(PROGRAMNAME) ../bin

cleandist:
	rm -f *.o filereader $(PROGRAMNAME) *~
