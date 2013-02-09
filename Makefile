HEADERS := graph_bgl_policy.hpp gspan.hpp 
SOURCES := main.cpp

CFLAGS := -p -O0 -g -Wall
CFLAGS += -DDEBUG_CHECK_GRAPH_LABEL #-DNDEBUG
CFLAGS += -DDEBUG_PRINT

all: closegraph

gspan: ${SOURCES} ${HEADERS}
	g++ -DGSPAN_FUNCTION=gspan      ${CFLAGS} ${SOURCES} -o gspan

closegraph: ${SOURCES} ${HEADERS}
	g++ -DGSPAN_FUNCTION=closegraph ${CFLAGS} ${SOURCES} -o closegraph

clean:
	rm -f *.o *.ii *.s gspan closegraph
