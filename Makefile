CC=g++ -O3 -Wall
INCLUDE= -I /home/yerui/src/alglib-3.14.0/src 
OBJS=main.o

twoPois: twoPois.cc $(OBJS)
	$(CC) $< -o $@ $(OBJS) /home/yerui/src/alglib-3.14.0/src/*.o $(INCLUDE)

main.o: main.cc main.h
	$(CC) -c $< -o $@ $(INCLUDE)

.PHONY: clean

clean:
	-rm twoPois *.o
