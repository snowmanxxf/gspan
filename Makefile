
SOURCES := main.cpp
HEADERS := gspan.hpp

gspan: $(SOURCES) $(HEADERS)
	g++ -O0 -g -Wall $(SOURCES) -o gspan