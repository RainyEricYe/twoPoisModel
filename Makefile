twoPois: main.cc main.h
	g++ main.cc -o twoPois /home/yerui/src/alglib-3.14.0/src/*.o -I /home/yerui/src/alglib-3.14.0/src

.PHONY: clean

clean:
	-rm twoPois
