CFLAGS = -Wall 
CC = gcc

objects = geneticPerm.o tour.o

genetic: $(objects)
	$(CC) $(CFLAGS) -o genetic $(objects)

geneticPerm.o: geneticPerm.c geneticPerm.h tour.h
tour.o: tour.c tour.h

.PHONY : clean
clean: 
	rm genetic $(objects)
