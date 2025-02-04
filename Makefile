CC     := mpicc
CFLAG  := -std=c99 -Wall -Wextra -O3 -DNDIMS=3
INC    := -Iinclude -ISimpleDecomp/include -ISimpleNpyIO/include
LIB    := -lfftw3 -lm
SRCDIR := src SimpleDecomp/src SimpleNpyIO/src
OBJDIR := obj
SRCS   := $(shell find $(SRCDIR) -type f -name *.c)
OBJS   := $(patsubst %.c,obj/%.o,$(SRCS))
DEPS   := $(patsubst %.c,obj/%.d,$(SRCS))
OUTDIR := output
TARGET := a.out

help:
	@echo "all     : create \"$(TARGET)\""
	@echo "clean   : remove \"$(TARGET)\" and object files under \"$(OBJDIR)\""
	@echo "output  : create \"$(OUTDIR)\" to store output"
	@echo "datadel : clean-up \"$(OUTDIR)\""
	@echo "help    : show this message"

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAG) -o $@ $^ $(LIB)

$(OBJDIR)/%.o: %.c
	@if [ ! -e $(dir $@) ]; then \
		mkdir -p $(dir $@); \
	fi
	$(CC) $(CFLAG) -MMD $(INC) -c $< -o $@

clean:
	$(RM) -r $(OBJDIR) $(TARGET)

output:
	@if [ ! -e $(OUTDIR)/log ]; then \
		mkdir -p $(OUTDIR)/log; \
	fi
	@if [ ! -e $(OUTDIR)/save ]; then \
		mkdir -p $(OUTDIR)/save; \
	fi
	@if [ ! -e $(OUTDIR)/stat ]; then \
		mkdir -p $(OUTDIR)/stat; \
	fi

datadel:
	$(RM) -r $(OUTDIR)/log/*
	$(RM) -r $(OUTDIR)/save/*
	$(RM) -r $(OUTDIR)/stat/*

-include $(DEPS)

.PHONY : all clean output datadel help

