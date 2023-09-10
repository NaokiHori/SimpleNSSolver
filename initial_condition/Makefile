help:
	@echo "output  : make directory to store NPY files"
	@echo "datadel : remove NPY files"
	@echo "help    : show this message"

output:
	@if [ ! -e output ]; then \
		mkdir output; \
	fi

datadel:
	$(RM) output/*.npy

.PHONY : help output datadel

