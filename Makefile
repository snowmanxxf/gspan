HEADERS := graph_policy.hpp gspan.hpp
SOURCES := main.cpp

MACRO   := -DDEBUG_CHECK_GRAPH_LABEL
#MACRO   += -DDEBUG_PRINT

all: gspan closegraph

gspan: ${SOURCES} ${HEADERS}
	g++ -DGSPAN_FUNCTION=gspan      ${MACRO} -O3 -Wall ${SOURCES} -o gspan

closegraph: ${SOURCES} ${HEADERS}
	g++ -DGSPAN_FUNCTION=closegraph ${MACRO} -O3 -Wall ${SOURCES} -o closegraph

clean:
	rm -f *.o *.ii *.s grviz/*
