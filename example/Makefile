CC=g++  -Wall
INCLUDE= -I /home/yerui/src/alglib-3.14.0/src  -I /home/yerui/src/SeqLib/ -I /home/yerui/src/SeqLib/htslib/ /home/yerui/src/SeqLib/bin/libseqlib.a /home/yerui/src/SeqLib/bin/libhts.a /home/yerui/src/SeqLib/bin/libbwa.a /home/yerui/src/SeqLib/bin/libfml.a -I /home/yerui/anaconda3/include  -L /home/yerui/anaconda3/lib
LIBS= -llzma -lbz2 -L. -lz -lpthread

OBJS=main.o

sens: sens.cc $(OBJS)
	$(CC) $< -o $@ $(OBJS) /home/yerui/src/alglib-3.14.0/src/*.o $(INCLUDE) $(LIBS)

main.o: ../main.cc ../main.h
	$(CC) -c $< -o $@ $(INCLUDE)

Option.o: Option.h
	$(CC) -c $< -o $@

.PHONY: clean

clean:
	-rm sens *.o
