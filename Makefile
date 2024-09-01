# clkcomb - Clock and phase bias products Combination
# Copyright (C) 2021 Yuanxin Pan
# All rights reserved.

DIR_GUARD = @mkdir -p $(@D)
BIN_PATH = ./bin/
OBJ_PATH = $(BIN_PATH)/obj/

SUBDIR = modules \
	     apps

.PHONY : all $(SUBDIR) clean

all : $(SUBDIR)

apps utils : modules

$(SUBDIR) :
	make -j4 -C src/$@

clean :
	-rm -rf $(BIN_PATH)
