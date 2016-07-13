include ../make.config

BUILDS = original_code tapenade #openad

build:${BUILDS}
	$(foreach d, $^, make -C $(d) build;)

clean:${BUILDS}
	$(foreach d, $^, make -C $(d) clean;)
