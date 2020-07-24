###############################################################################
# main/Makefile                                                               #
#                                                                             #
# Makefile for the main folder of ader program. This is treated special and	#
# does not include the Makefile.subdirs. Here the ader executable is linked.  #
###############################################################################

NAME := main

OBJS := 
OBJS += main.o 
OBJS += ader_dg.o
OBJS += system.o
OBJS += func_pointer.o
OBJS += iteration.o
OBJS += rk_predictor.o


LIBS         = $(wildcard $(LIBDIR)/*.a)
OBJDIR       = $(LIBDIR)/obj/$(NAME)
OBJSWITHPATH = $(addprefix $(OBJDIR)/,$(OBJS))



$(EXEDIR)/$(EXE): $(OBJSWITHPATH) $(LIBS)
	$(CC) $(CFLAGS) -o $@ $(OBJSWITHPATH) $(libsys) $(LIBS) $(LIBS) $(LIBS) -L$(LIBDIR) $(libsys)

$(OBJDIR)/%.o: $(NAME).h %.c Makefile
	$(CC) $(CFLAGS) $(incsglib) -o $@ -c $*.c 

