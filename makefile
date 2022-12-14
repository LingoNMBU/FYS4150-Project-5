all:	compile1 run1 compile2 run2

compile1:
	g++ -O2 prob7.cpp src/wavebox.cpp -I include -o prob7.exe -larmadillo
run1:
	./prob7.exe
compile2:
	g++ -O2 prob8.cpp src/wavebox.cpp -I include -o prob8.exe -larmadillo
run2:
	./prob8.exe