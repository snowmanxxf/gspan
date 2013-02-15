HEADERS := graph_bgl_csr_policy.hpp graph_bgl_adjl_policy.hpp gspan.hpp
SOURCES := main.cpp

CFLAGS := -O3 -g -p -Wall
CFLAGS += -DGRAPH_ADJL
CFLAGS += -DDEBUG_CHECK_GRAPH_LABEL #-DNDEBUG
#CFLAGS += -DDEBUG_PRINT

all: closegraph_st closegraph_mt

closegraph_st: ${SOURCES} ${HEADERS} closegraph_st.hpp
	g++ -DCLOSEGRAPH_ST ${CFLAGS} ${SOURCES} -o closegraph_st

closegraph_mt: ${SOURCES} ${HEADERS}  closegraph_mt.hpp
	g++ -DCLOSEGRAPH_MT -lboost_thread ${CFLAGS} ${SOURCES} -o closegraph_mt

#gspan: ${SOURCES} ${HEADERS}
#	g++ ${CFLAGS} ${SOURCES} -o gspan

clean:
	rm -f *.o *.ii *.s gspan closegraph_mt closegraph_st
