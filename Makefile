
CFLAGS := -O3 -g -Wall -march=amdfam10 -save-temps #-fno-inline

CFLAGS += -falign-loops -fprefetch-loop-arrays -freg-struct-return
CFLAGS += -fschedule-insns -fsched-pressure

#CFLAGS +=  -fprofile-generate
#CFLAGS +=  -fprofile-use

CFLAGS += -DDEBUG_CHECK_GRAPH_LABEL 
CFLAGS += -DGSPAN_WITH_STATISTICS
CFLAGS += -DGSPAN_TRACE
#CFLAGS += -DCHECK_MODE
#CFLAGS += -DWITH_CHECKS
#CFLAGS += -DNDEBUG
#CFLAGS += -DTYPE_CHECK
CFLAGS += -DUSE_ASM

# test programm
closegraph: main.cpp misc.hpp read_input.cpp gspan.hpp libgspan.a
	g++ ${CFLAGS} main.cpp read_input.cpp -L. -lgspan -o closegraph


# gspan static library
libgspan.a: gspan.o gspan_graph.o gspan_allocator.o
	ar rcs libgspan.a gspan.o gspan_graph.o gspan_allocator.o

gspan.o: gspan.cpp gspan.hpp gspan_graph.o gspan_allocator.o
	g++ ${CFLAGS} -c gspan.cpp -o gspan.o

gspan_graph.o: gspan_graph.hpp gspan_graph.cpp
	g++ $(CFLAGS) -c gspan_graph.cpp -o gspan_graph.o

gspan_allocator.o: gspan_allocator.hpp gspan_allocator.cpp
	g++ ${CFLAGS} -c gspan_allocator.cpp -o gspan_allocator.o

VALGRIND_FILES := cachegrind.out.* callgrind.out.* massif.out.*

clean:
	rm -f *.o *.ii *.s *.gcda gmon.out libgspan.a closegraph ${VALGRIND_FILES}
