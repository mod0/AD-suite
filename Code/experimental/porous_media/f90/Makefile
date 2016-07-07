include ../make.config

BUILDS = original_code openad tapenade

build:${BUILDS}
	$(foreach d, $^, make -c $(d) build;)

clean:${BUILDS}
	$(foreach d, $^, make -c $(d) clean;)
