PROG := inpaint
SRCS := $(wildcard *.c) minitrace/minitrace.c
OBJS := $(SRCS:%.c=%.o)
DEPS := $(SRCS:%.c=%.d)

CC := gcc
CCFLAGS := -D_DEFAULT_SOURCE -std=c11 -O0 -g -Wall -DMTR_ENABLED
INCLUDEPATH := -I/usr/local/include 
LIBPATH := -L/usr/local/lib 
LIBS := 
STATICS := 
LNFLAGS := 


all: $(DEPENDS) $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CCFLAGS) -o $@ $^ $(LIBPATH) $(STATICS) $(LIBS) $(LNFLAGS)

.c.o:
	$(CC) $(CCFLAGS) $(INCLUDEPATH) -MMD -MP -MF $(<:%.c=%.d) -c $< -o $(<:%.c=%.o)


.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) $(DEPS)

-include $(DEPS)
