
SOURCES := main.cpp
HEADERS := gspan.hpp graph_ops.hpp

gspan: $(SOURCES) $(HEADERS)
	g++ -DDEBUG_PRINT -O2 -p -Wall $(SOURCES) -o gspan
