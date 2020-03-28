CC = g++ 
CFLAGS = -O4 -Wno-deprecated -Wall

HDRS = eclat.h timetrack.h calcdb.h eqclass.h stats.h\
	maximal.h chashtable.h
OBJS = calcdb.o eqclass.o stats.o maximal.o eclat.o enumerate.o\
	chashtable.o
LIBS = 
TARGET = schism

default: $(TARGET)

clean:
	rm -rf *~ *.o $(TARGET) dataGen convertData

schism: $(OBJS) $(HDRS)
	$(CC) $(CFLAGS) -o eclat $(OBJS) $(LIBS)

dna_eclat: $(OBJS) $(HDRS)
	$(CC) $(CFLAGS) -DDNA -o dna_eclat $(OBJS) $(LIBS)

makebin: 
	$(CC) $(CFLAGS) -o makebin makebin.cc

dataGen: 
	$(CC) $(CFLAGS) -o dataGen dataGen.cc

convertData: 
	$(CC) $(CFLAGS) -o convertData convertData.cc

.SUFFIXES: .o .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $<


# dependencies
# $(OBJS): $(HDRS)
eclat.o: $(HDRS)
enumerate.o: $(HDRS)
calcdb.o: calcdb.h eclat.h eqclass.h
eqclass.o: eqclass.h eclat.h calcdb.h
maximal.o: maximal.h calcdb.h eqclass.h
stats.o: stats.h
chashtable.o: chashtable.h eclat.h
