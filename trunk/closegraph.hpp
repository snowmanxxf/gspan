#ifndef CLOSEGRAPH_H_
#define CLOSEGRAPH_H_

#include "gspan.hpp"

namespace graph_alg
{
    namespace CloseGraph
    {
	// *****************************************************************************
	//                          functions
	// *****************************************************************************

	template<class TGraphIterator, class GraphOps, class Output>
	void closegraph(TGraphIterator tg_begin, TGraphIterator tg_end, int minsup,
			const GraphOps& ops, Output& result)
	{
	}
    }
}

#endif
