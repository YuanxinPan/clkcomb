
CLK_OBJ := main.o \
           utils.o \
           config.o \
           AnalyseCenter.o \

OBJ += $(CLK_OBJ)
CLK_OBJ := $(addprefix $(OBJ_PATH), $(CLK_OBJ))

$(BIN_PATH)clkcomb : $(CLK_OBJ) $(LIB)
	$(CC) $(CPPFALGS) $(DEBUG) $(CLK_OBJ) $(LIB) -o $(BIN_PATH)clkcomb

