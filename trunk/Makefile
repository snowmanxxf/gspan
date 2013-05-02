HEADERS := graph_bgl_csr_policy.hpp graph_bgl_adjl_policy.hpp gspan.hpp
SOURCES := main.cpp

CFLAGS := -O3 -g -Wall -march=amdfam10
CFLAGS += -DGRAPH_ADJL
CFLAGS += -DDEBUG_CHECK_GRAPH_LABEL 
CFLAGS += -DGSPAN_WITH_STATISTICS
CFLAGS += -DGSPAN_TRACE
#CFLAGS += -DCHECK_MODE
CFLAGS += -DDEBUG_PRINT
#CFLAGS += -DNDEBUG

closegraph: main2.cpp gspan2.hpp gspan2.cpp misc.hpp read_input.cpp gspan_graph.hpp
	g++ ${CFLAGS} -DGSPAN_FUNCTION=closegraph main2.cpp gspan2.cpp read_input.cpp -o closegraph


clean:
	rm -f *.o *.ii *.s gspan closegraph_mt closegraph_st

