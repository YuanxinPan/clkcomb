
ATT_OBJ := att.o

OBJ += $(ATT_OBJ)
ATT_OBJ := $(addprefix $(OBJ_PATH), $(ATT_OBJ))

$(BIN_PATH)att : $(ATT_OBJ) $(LIB)
	$(CC) $(CPPFALGS) $(DEBUG) $(ATT_OBJ) $(LIB) -o $(BIN_PATH)att

