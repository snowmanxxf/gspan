#ifndef READ_INPUT_H_
#define READ_INPUT_H_

#include <vector>
#include <map>
#include <list>
#include <string>
#include <iostream>

struct InputGraph
{
    std::string name;
    std::map<int, std::string> vl;

    struct E
    {
	int from, to;
	std::string el;
    };

    std::vector<E> edges;
};

void read_input(std::list<InputGraph>& igl, std::istream& is);


#endif
