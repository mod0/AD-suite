include make.config

define HELP_BODY
This makefile has the following Options and
corresponding values that can be passed to 
it from the command line while calling make.

+-------------------------------------+
|              OPTIONS                |
+-------------------------------------+
LANG     =   <ALL/F77/>
DERIV    =   <ALL/TAPENADE/OPENAD>
endef

export HELP_BODY

# Space separated list of folders to build
BUILDS = f77

.PHONY: all build help clean

all: build

build: ${BUILDS}
	$(foreach d, $^, make -C $(d) build;)

help:
	@echo "$$HELP_BODY"

clean: ${BUILDS}
	$(foreach d, $^, make -C $(d) clean;)

