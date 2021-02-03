DIR_GUARD = @mkdir -p $(@D)
BIN_PATH = ../bin/
OBJ_PATH = ../bin/obj/utils/

EXE := fip2snx grg2snx sgg2snx zqy2snx
TARGET := $(addprefix $(BIN_PATH), $(EXE))
all : $(TARGET)

CC = g++
CPPFALGS = -Wall -g -O3 -I../include -std=c++11
DEBUG = 

SRC := $(wildcard *.cpp)
OBJ := $(SRC:.f90=.o)
OBJ := $(addprefix $(OBJ_PATH), $(OBJ))

$(OBJ_PATH)%.o : %.cpp
	$(DIR_GUARD)
	$(CC) $(CPPFALGS) $(DEBUG) -c $< -o $@
	@$(CC) -MM $(CPPFALGS) $< | \
    sed 's,\($*\)\.o[ :]*,$@ : ,g' > $(@:.o=.d)

$(BIN_PATH)fip2snx : $(OBJ_PATH)fip2snx.o
	$(CC) $(CPPFALGS) $(DEBUG) $< -o $@

$(BIN_PATH)grg2snx : $(OBJ_PATH)grg2snx.o
	$(CC) $(CPPFALGS) $(DEBUG) $< -o $@

$(BIN_PATH)sgg2snx : $(OBJ_PATH)sgg2snx.o
	$(CC) $(CPPFALGS) $(DEBUG) $< -o $@

$(BIN_PATH)zqy2snx : $(OBJ_PATH)zqy2snx.o
	$(CC) $(CPPFALGS) $(DEBUG) $< -o $@

.PHONY : clean
clean :
	-rm $(TARGET) $(OBJ_PATH)*

-include $(OBJ:.o=.d)

