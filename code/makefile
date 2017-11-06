ROOT_FLAGS = `root-config --cflags`
ROOT_LIBS  = `root-config --libs`
CFLAGS		 = -I.
RCXX = $(CFLAGS) $(ROOT_FLAGS) $(ROOT_LIBS) $(LIBS)

INC_PATH = -I/usr/include -Iinc
LIB_PATH = -L/usr/lib -L/usr/local/lib 
LIBS     = -lboost_system -lMinuit -lpthread
           

# compiler flags:
#    -g    adds debugging information to the executable file
#    -Wall turns on most, but not all, compiler warnings

DBG    = -g -O3 -Wall 

CC     = g++ $(DBG) $(LIB_PATH) $(LIBS) $(INC_PATH) -lpthread -std=c++11
MIN_CC = g++ $(DBG) $(INC_PATH) $(CFLAGS) $(ROOT_FLAGS) -std=c++11

BIN_DIR = bin
OBJ_DIR = obj
SRC_DIR = src
INC_DIR = inc

all: FastRotationToyModel FastRotationToyModelCalo FFT fillDigitizeFRS fillWshape 

fillWshape:	$(SRC_DIR)/fillWshape.cpp
	$(CC) -o $(BIN_DIR)/fillWshape $(SRC_DIR)/fillWshape.cpp $(RCXX)

fillDigitizeFRS:	$(SRC_DIR)/fillDigitizeFRS.cpp
	$(CC) -o $(BIN_DIR)/fillDigitizeFRS $(SRC_DIR)/fillDigitizeFRS.cpp $(RCXX)

FastRotationToyModel: $(SRC_DIR)/FastRotationToyModel.cpp
	$(CC) -o $(BIN_DIR)/FastRotationToyModel $(SRC_DIR)/FastRotationToyModel.cpp $(RCXX)

FastRotationToyModelCalo: $(SRC_DIR)/FastRotationToyModelCalo.cpp
	$(CC) -o $(BIN_DIR)/FastRotationToyModelCalo $(SRC_DIR)/FastRotationToyModelCalo.cpp $(RCXX)

FFT: $(SRC_DIR)/FFT.cpp
	$(CC) -o $(BIN_DIR)/FFT $(SRC_DIR)/FFT.cpp $(RCXX)

clean:
	rm -f $(OBJ_DIR)/*.o $(BIN_DIR)/*
