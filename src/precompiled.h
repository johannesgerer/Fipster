#pragma once


//for assert.h / cassert
#ifndef DEBUG
// #define NDEBUG 1
#endif
 
#include "src/precompiled_auto.h"

#ifdef __INTEL_COMPILER
#define nullptr NULL
#endif

// in this file include headers, that do no change much during
// development. they will be included twice but usually contain
// '#pragma once'. headers that should alwyas be precompiled are
// placed within autocompilation regions and then extracted by the
// make file (this requires a make cleanprecomp).

/** Select the sub grid for the true time-stepping variables. 
	It has to be larger that 0 to have one at least a one point thick boundary
*/
#define J_SGS 1 

#if 1

#include "integer_types.h"

#include "named_counter.h"
#include "auto_connecting_body.h"
// #include "memoized_node_factory.h"
#include "rooted_graph.h"

/*
#include "spacetime_field.h"

#include "grid.h"
#include "grid_fields.h"
#include "grid_indices.h"
#include "grid_iterators.h"
*/

#include "logging.h" 


#endif
