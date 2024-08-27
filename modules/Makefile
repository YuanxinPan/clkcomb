# clkcomb - Clock and phase bias products Combination
# Copyright (C) 2021 Yuanxin Pan
# All rights reserved.

DIR_GUARD = @mkdir -p $(@D)
BIN_PATH = ../bin/
OBJ_PATH = ../bin/obj/modules/

TARGET := $(BIN_PATH)libmods.a
all : $(TARGET)

CC = g++
CPPFALGS = -Wall -g -O3 -I../include -std=c++11
DEBUG = 

MODULE := chrono \
		  coord \
		  io \
		  rinex
VPATH := $(MODULE)
OBJ := 
include $(addsuffix /Module.mk, $(MODULE))
OBJ :=  $(addprefix $(OBJ_PATH), $(OBJ))

$(TARGET) : $(OBJ)
	$(DIR_GUARD)
	ar -rv $(TARGET) $(OBJ)

$(OBJ_PATH)%.o : %.cpp
	$(DIR_GUARD)
	$(CC) $(CPPFALGS) $(DEBUG) -c $< -o $@
	@$(CC) -MM $(CPPFALGS) $< | \
    sed 's,\($*\)\.o[ :]*,$@ : ,g' > $(@:.o=.d)

.PHONY : clean
clean :
	-rm $(TARGET) $(OBJ_PATH)*

-include $(OBJ:.o=.d)

