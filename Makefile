DIR_GUARD = @mkdir -p $(@D)
BIN_PATH = ./bin/
OBJ_PATH = ./bin/obj/

SUBDIR = modules \
	     apps

.PHONY : all $(SUBDIR) clean

all : $(SUBDIR)

apps utils : modules

$(SUBDIR) :
	make -j4 -C $@

clean :
	-rm -r $(BIN_PATH)
