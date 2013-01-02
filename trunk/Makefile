HEADERS := graph_policy.hpp gspan.hpp edge_iterator.hpp
SOURCES := main.cpp

CFLAGS := -O3 -p -g -Wall
CFLAGS += -DDEBUG_CHECK_GRAPH_LABEL
#CFLAGS += -DDEBUG_PRINT

all: gspan closegraph

gspan: ${SOURCES} ${HEADERS}
	g++ -DGSPAN_FUNCTION=gspan      ${CFLAGS} ${SOURCES} -o gspan

closegraph: ${SOURCES} ${HEADERS}
	g++ -DGSPAN_FUNCTION=closegraph ${CFLAGS} ${SOURCES} -o closegraph

clean:
	rm -f *.o *.ii *.s gspan closegraph
