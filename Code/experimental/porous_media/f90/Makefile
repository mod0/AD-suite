include ../make.config

BUILDS = original_code #openad tapenade

build:${BUILDS}
	$(foreach d, $^, make -C $(d) build;)

clean:${BUILDS}
	$(foreach d, $^, make -C $(d) clean;)
	rm -i *.png *.jpg *.txt
