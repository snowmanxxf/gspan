
SOURCES := main.cpp
HEADERS := gspan.hpp graph_ops.hpp

gspan: $(SOURCES) $(HEADERS)
	g++ -O2 -Wall $(SOURCES) -o gspan
