#include "precompiled.h"

#include "named_counter.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <pretty_printer.h>
#endif


namespace fipster { 


	void named_counter::print()
	{
		print_line(cout,get_singleton().counters);
	}

}
