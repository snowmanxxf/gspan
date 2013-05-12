HEADERS := graph_bgl_csr_policy.hpp graph_bgl_adjl_policy.hpp gspan.hpp
SOURCES := main.cpp

CFLAGS := -O3 -g -Wall -march=amdfam10 -save-temps
CFLAGS += -DGRAPH_ADJL
CFLAGS += -DDEBUG_CHECK_GRAPH_LABEL 
CFLAGS += -DGSPAN_WITH_STATISTICS
CFLAGS += -DGSPAN_TRACE
#CFLAGS += -DCHECK_MODE
CFLAGS += -DDEBUG_PRINT
CFLAGS += -DNDEBUG
#CFLAGS += -DTYPE_CHECK
CFLAGS += -DUSE_ASM


closegraph: main.cpp misc.hpp read_input.cpp gspan.o gspan_allocator.o
	g++ ${CFLAGS} main.cpp read_input.cpp gspan.o gspan_allocator.o -o closegraph

gspan.o: gspan.cpp gspan.hpp gspan_graph.hpp gspan_allocator.o
	g++ ${CFLAGS} -c gspan.cpp -o gspan.o

gspan_allocator.o: gspan_allocator.hpp gspan_allocator.cpp
	g++ ${CFLAGS} -c gspan_allocator.cpp -o gspan_allocator.o

clean:
	rm -f *.o *.ii *.s closegraph
