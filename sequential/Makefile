CC = gcc
LD = gcc
CFLAGS = -Wall
LDFLAGS = -lm
RM = /bin/rm -f
OBJS = lu_seq.o
EXECUTABLE = lu_seq

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(OBJS) -o $(EXECUTABLE) $(LDFLAGS)


$(OBJS):  lu_seq.c
	$(CC) $(CFLAGS) -c lu_seq.c 


clean:
	$(RM) $(EXECUTABLE) $(OBJS)