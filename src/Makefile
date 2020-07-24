###############################################################################
# @file Makefile                                                              #
# @author mapi 03/15											                           #
# @brief build ader_dg program		                                          #
###############################################################################

###############################################################################
# Include host configuration                                                  #
###############################################################################

include hosts.mk

###############################################################################
# Set variables                                                               #
###############################################################################

include config.mk

###############################################################################
# Set libraries to compile                                                     #
###############################################################################

include libs.mk

###############################################################################
# Export variables to make them visisble to the submake files                 #
###############################################################################

export

###############################################################################
# Targets                                                                     #
###############################################################################

# Main target
# Process sub directories first before going to main and building bamps
all:
	@for X in $(libs); do mkdir -p $(OBJDIR)/$$X; done
	@for X in $(libs); do $(MAKE)  -C $$X; done
	@mkdir -p $(OBJDIR)/$(MAIN)
	$(MAKE)  -C $(MAIN)

# init
init:
	gcc -g $(MSC)/initial_data/averages.c $(MSC)/initial_data/TOV.c -o $(MSC)/initial_data/run_average -lm
	-rm -r $(MATHDIR)/initial_data/*

# conv
conv:
	gcc -C $(CONVDIR)/convergence.c -o $(CONVDIR)/run_TEST -lm
	gcc -C $(CONVDIR)/newton_ana.c -o $(CONVDIR)/run_burger_ana -lm
	gcc -C $(CONVDIR)/conv_evolve.c -o $(CONVDIR)/run_con_TOV_evolve -lm
	gcc -C $(CONVDIR)/conv_evolve_cont.c -o $(CONVDIR)/run_con_TOV_evolve_cont -lm
	gcc -C $(CONVDIR)/conv_origin.c -o $(CONVDIR)/run_con_TOV_origin -lm

# clean
clean: 
	-rm -r $(LIBDIR)/*
	-rm -r $(EXEDIR)/$(EXE)

# new
new: clean all init conv

