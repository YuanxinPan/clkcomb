DIR_GUARD = @mkdir -p $(@D)
BIN_PATH = ../bin/
OBJ_PATH = ../bin/obj/apps/

MODULE := clkcomb att
          
TARGET := $(addprefix $(BIN_PATH), $(MODULE))
all : $(TARGET)

CC = g++
CPPFALGS = -Wall -g -O3 -I../include/ -std=c++11
DEBUG = 
LIB = $(BIN_PATH)libmods.a

VPATH := $(MODULE)
OBJ :=
include $(addsuffix /Module.mk, $(MODULE))
OBJ := $(addprefix $(OBJ_PATH), $(OBJ))

$(OBJ_PATH)%.o : $(SRC_PATH)%.cpp
	$(DIR_GUARD)
	$(CC) $(CPPFALGS) $(DEBUG) -c $< -o $@
	@$(CC) -MM $(CPPFALGS) $< | \
    sed 's,\($*\)\.o[ :]*,$@ : ,g' > $(@:.o=.d)

.PHONY : clean
clean :
	-rm $(TARGET) $(OBJ_PATH)*

-include $(OBJ:.o=.d)

