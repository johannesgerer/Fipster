#include "precompiled.h"
#include "finite_difference_weights_config.h"
#include "PDE.h"
#include "configuration_parser.h"

namespace fipster { namespace finite_difference_weights {



	// #################################################
	// ##############   PHASE 1:			############
	// ##############   CONFIGURATION       ############
	// #################################################

	config_t::config_t( const ptree& pt,objs_t objs )
		: id(betterGet<string>(pt,"<xmlattr>.id","finite difference weights"))
	{

		string stype = pt.get<string>("stencilType");

		if(stype == "Full") stencil_type = StencilFull;
		else if(stype == "NoDiagonals") stencil_type = StencilNoDiagonals;
		else if(stype == "OneDiagonal") stencil_type = StencilOneDiagonal;
		else if(stype == "TwoDiagonals") stencil_type = StencilTwoDiagonals;
		else BOOST_THROW_EXCEPTION(runtime_error(stype+" is no valid Stencil Type"));

		partial_derivative_order = pt.get<int>("partialDerivativeOrder");

		if(partial_derivative_order>2)
			higher_order_weight = pt.get<double>("higherOrderWeight");

		auto& pt_pde=betterChild(pt,"PDE",id);

		stype = pt_pde.get<string>("<xmlattr>.type");
		if(stype == "BSIso"){
			PDE.reset(new BSiso(id,pt_pde,stype));
		}else if(stype == "BSIsoGeom"){
			PDE.reset(new BSiso_geom_average(id,pt_pde,objs));
		}else BOOST_THROW_EXCEPTION(runtime_error
																(stype+" is no valid PDE Type"));
	}

	// #################################################
	// ##############   PHASE 2:			############
	// ##############   Node Factory        ############
	// #################################################
	

	



}}
