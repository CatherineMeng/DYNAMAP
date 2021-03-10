#
# generate dependencies automatically
#

# 
CC = gcc
CCFLAGS = -g -Wall -O2
LINK.o = gcc -g -O2 

# 
libfiles = interface.c bf_pbqp.c h_pbqp.c vec_mat.c  
csources = $(libfiles) pbqp_solve.c pbqp_validate.c
exe = pbqp_solve pbqp_validate 

# phony targets
all: $(exe)

clean: 
	rm -f *.o *.d $(exe)
install: all
	cp $(exe) ../bin


# dependence rules
%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

%.d: %.c
	@echo "create dependence for $<";$(SHELL) -ec '$(CC) -MG -MM $(CCFLAGS) $< \
	| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
	[ -s $@ ] || rm -f $@'

include $(csources:.c=.d)

#############################################################################
# targets 

pbqp_solve: $(libfiles:.c=.o)  pbqp_solve.o

pbqp_validate: $(libfiles:.c=.o)  pbqp_validate.o

