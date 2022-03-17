CC	= gcc
CFLAGS  = -O3 -Wall
LIBS    = -lm
VPATH   = mican/mican_2015.03.09/src
SRC     = $(shell ls $(VPATH)/*.c)
HEADERS = $(shell ls $(VPATH)/*.h)
OBJS    = $(SRC:.c=.o)

TARGET  = mican_2015.03.09

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

$(OBJS): $(HEADERS)

clean:
	rm $(OBJS) $(TARGET) $(TARGETLIB)
