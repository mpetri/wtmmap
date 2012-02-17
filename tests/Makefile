INCLUDES	:= -I ./CppUnitLite -I ../include
COMMON		:= ./CppUnitLite/*.cpp test-main.cpp

all: clean rankbvTest wtTest run

rankbvTest:
	g++ -Wall -g -o rankbvTest $(INCLUDES) $(COMMON) ../src/rankbv.c rankbvTest.cpp

wtTest:
	g++ -Wall -g -o wtTest $(INCLUDES) $(COMMON) ../src/cbheap.c ../src/rankbv.c ../src/wt.c wtTest.cpp

run:
	./rankbvTest
	./wtTest

clean:
	rm -f ./rankbvTest
	rm -f ./wtTest 
