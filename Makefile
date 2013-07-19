
CXXFLAGS := -O3 -fPIC -g -Wall

CXXFLAGS += -falign-loops -fprefetch-loop-arrays -freg-struct-return
CXXFLAGS += -fschedule-insns -fsched-pressure

#CXXFLAGS +=  -fprofile-generate
#CXXFLAGS +=  -fprofile-use

CXXFLAGS += -DDEBUG_CHECK_GRAPH_LABEL 
CXXFLAGS += -DGSPAN_WITH_STATISTICS
CXXFLAGS += -DGSPAN_TRACE
#CXXFLAGS += -DCHECK_MODE
#CXXFLAGS += -DWITH_CHECKS
#CXXFLAGS += -DNDEBUG
#CXXFLAGS += -DTYPE_CHECK
CXXFLAGS += -DUSE_ASM
#CXXFLAGS += -DET_TEST_PATH_MIN


# test programms
closegraph_lemon: main_lemon.cpp gspan_lemon.hpp intlabelmap.hpp libgspan.a
	g++ ${CXXFLAGS} main_lemon.cpp -L. -lgspan -lemon -o closegraph_lemon

closegraph: main.cpp read_input.cpp misc.hpp misc.cpp gspan.hpp libgspan.a
	g++ ${CXXFLAGS} main.cpp read_input.cpp misc.cpp -L. -lgspan -o closegraph


# gspan static library
libgspan.a: gspan.o gspan_graph.o gspan_allocator.o
	ar rcs libgspan.a gspan.o gspan_graph.o gspan_allocator.o

gspan.o: gspan.cpp gspan.hpp gspan_graph.o gspan_allocator.o
	g++ ${CXXFLAGS} -c gspan.cpp -o gspan.o

gspan_graph.o: gspan_graph.hpp gspan_graph.cpp
	g++ $(CXXFLAGS) -c gspan_graph.cpp -o gspan_graph.o

gspan_allocator.o: gspan_allocator.hpp gspan_allocator.cpp
	g++ ${CXXFLAGS} -c gspan_allocator.cpp -o gspan_allocator.o

VALGRIND_FILES := cachegrind.out.* callgrind.out.* massif.out.*

clean:
	rm -f *.o *.ii *.s *.gcda gmon.out libgspan.a closegraph ${VALGRIND_FILES}
