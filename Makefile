INCLUDES	:= -I ./include 

all: clean tests

tests:
	make -C tests

clean:
	make clean -C tests

.PHONY: tests
