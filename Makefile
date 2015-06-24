all: TreeReader.run

TreeReader.run: TreeReader.o 
	g++ -o TreeReader.run TreeReader.o  `root-config --libs`

TreeReader.o: TreeReader.C
	g++ -I`root-config --incdir` -c -g TreeReader.C

clean:
	rm TreeReader.run TreeReader.o
