CC = gcc
CFLAGS = -Wall -fopenmp -O1 -g

TARGET = app
OUTPUT_FILE = spikes.txt
VISUALIZATION_FILE = spikes.npy

OBJ_DIR = ./obj
SRC_DIR = ./src

SRCS := $(shell find $(SRC_DIR) -name '*.c')
OBJS := $(SRCS:%=$(OBJ_DIR)/%.o)

LIBS = -lm

RUN_FLAGS = -r


$(TARGET): clean $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS)

$(OBJ_DIR)/%.c.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

run:
	./$(TARGET) $(RUN_FLAGS)

.PHONY: clean

clean:
	rm -rf $(OBJ_DIR)/* $(TARGET) $(OUTPUT_FILE) $(VISUALIZATION_FILE)