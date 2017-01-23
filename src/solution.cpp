#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <fstream>
#include <array>
#endif

#include "solution.h"
#include "auto_connecting_body.h"
#include "expectation_value_node_factory.h"
#include "rooted_graph.h"
#include "global_config.h"
#include "Math.h"
#include "option.h"
#include "grid.h"
#include "spacetime_field.h"
#include "runtime_join_node.h"
#include "task_queue.h"
#include "interpolation_hedging.h"

namespace fipster { namespace _solution {  

	struct result_dumper_body : 
			auto_connecting_body<result_dumper_body //body_t
								,int				//result_t
								,vector<spacetime_field<J_SGS>::result_sptr_t> //future_t
			>
	{		
		solution& sol;
		btime time;
		result_dumper_body(solution* sol,btime time)
			:sol(*sol),time(time)
		{}

		result_sptr_t operator()(const input_t& field_results){
			sol.not_done.erase(time);
			//logging 
			if(0)thread_logger()<<"OUTPUT ("<<sol<<","<<time<<") | START"<<endl;
      if(1)thread_logger()<<"OUTPUT ("<<sol.id<<" at time "<<
						 time<<" on grid "<<sol.grid_name<<") | START"<<endl;
			
			auto name = global_config::get().output_path+sol.id+"-"+sol.grid_name+"-"+toS(time);
			ofstream out (name+".csv");
			
			if(!out.good())
				FIPSTER_THROW_EXCEPTION(runtime_error("Failed to open "+name+".csv"));
			
			ofstream out2;
			if(sol.single_values.size()>0){
				out2.open(name+"-singleValues.csv");
				if(!out2.good())
				FIPSTER_THROW_EXCEPTION(runtime_error("Failed to open "+name+"-singleValues.csv"));
			}

			//grid_field<grid_index<J_SGS>,vector<string>> lines;

			const uint D=sol.grid->D;

			{
				auto field=field_results->begin();
				FIPSTER_ASSERT(field_results->size() == sol.fields.size());
				for(auto& f : sol.fields){
					// \todo CodeRefBadHack
					if(!f.hedging || !get<2>(*f.hedging))
						out<<"#"<<f.id<<": "<<(*field++)->meta_info<<endl;
				}
			}

			out<<sol<<endl;

			//output single values
			if(!sol.single_values.empty()){
				out2<<sol<<endl;			
				for(auto k=begin(sol.single_values);k!=sol.single_values.end();k++){
					FIPSTER_ASSERT(k->size()==sol.grid->D);
					if(k->size() > 1)
						FIPSTER_THROW_EXCEPTION(runtime_error("only 1D interpolation implemented!"));
					double &v = (*k)[0];
					out2<<v<<"\t";
					//output one field's value at a time:
					for(auto& j  : *field_results){
						out2<<boost::lexical_cast<string>(
							interpolate::oneD(
								sol.grid->coords[0].begin()+J_SGS,
								sol.grid->coords[0].end()-J_SGS,
								**j,0,v)
						);
						out2<<"\t";
					}
					out2<<"\n";
				}
			}

			auto gend=sol.grid->end<J_SGS>();
			grid_iterator<1,true> gbegin(*sol.grid);

			
			
			
			vector<boost::optional<options::bands_field> > bands2(sol.fields.size());
			vector<array<discretized_space_field<1>,4> > band_results(sol.fields.size());

			
			//interpolation of hedging positions
			vector<boost::optional<discretized_space_field<1> > > interpolated(sol.fields.size());

			const auto field=sol.fields.begin();
			const auto result=field_results->begin();

			for(uint j=0; j<sol.fields.size(); ++j)
				if(field[j].hedging){
					auto& h  = *(field[j].hedging);
					auto hg=get<0>(h);
						
					if(!get<2>(*field[j].hedging)){
						auto f = dynamic_pointer_cast<const option_field_t>(result[j]);
						if(f->additional.at(gbegin.index()).second){
							auto &i = *(interpolated[j] = discretized_space_field<1>(""));
							options::interpolate_hedging(options::interpolate_option(*f)
																					 , *sol.grid, *hg, gbegin, gend, i);
						}
					}else{
						//calculate hedging bands
						auto &bands3 = *(bands2[j] = options::bands_field("hallo"));
						auto &band_result = band_results[j];
						bands3.resize(*sol.grid);
						grid_iterator<0,true> hit(*hg);
						auto hend = hit.ind.end(*hg);
						for(;hit != hend; ++hit, ++j){
							auto f = dynamic_pointer_cast<const option_field_t>(result[j]);
							grid_iterator<J_SGS,true> state_it(*sol.grid);
							for(;state_it!=gend;++state_it){
								auto &bands = bands3.at(state_it.ind);
								auto optimal = f->additional.at(state_it.ind).second;
								if(!optimal) continue;
								if(*optimal < hit.ind)
									bands[1] = optimal;
								else if(*optimal > hit.ind)
									bands[0] = optimal;
								else{
									auto compare=[&](const grid_index<1>& c){
										auto c_op =  f->additional.at(c).second;
										if(c_op){
											if(*c_op < hit.ind){
												auto& b = bands[3];
												if(!b || *b < hit.ind) b = hit.ind;
											}else if(*c_op > hit.ind){
												auto& b = bands[2];
												if(!b || *b > hit.ind) b = hit.ind;
											}
										}
									};
									if(state_it.ind+1 != gend)   compare(state_it.ind + 1);
									if(state_it.ind != gbegin.index()) compare(state_it.ind - 1);
								}
							}
						}
						--j;
						for(uint i=0; i<2; i++){
							// break;
							options::interpolate_hedging(options::interpolate_band(bands3,i)
																					 , *sol.grid, *hg, gbegin, gend, band_result[i]);
						}
					}
				}
			
			//output states:
			grid_iterator<J_SGS,true> state_it(*sol.grid);
			for(;state_it!=gend;++state_it){
				for(uint j=0;j<D;j++)
					out<<boost::lexical_cast<string>(state_it.state()[j])<<"\t";				
				//output one field's value at a time:
				for(uint j=0; j<sol.fields.size(); ++j){
					auto grid1D=field[j].grid;
					if(!field[j].hedging || !get<2>(*field[j].hedging)){
						if(!grid1D)
							out<<boost::lexical_cast<string>(result[j]->at(state_it.ind));
						else
							out<<boost::lexical_cast<string>
								(interpolate::oneD
								 (grid1D->coords[0].begin()+J_SGS,
									grid1D->coords[0].end()-J_SGS,
									**result[j],0,
									geometric_average(state_it.state(),D))
								 );
						out<<"\t";
						auto f = dynamic_pointer_cast<const option_field_t>(result[j]);
						if(f){
							if(field[j].hedging){
								string res;
								if(0){ //use interpolation
									auto& h  = *(field[j].hedging);
									auto hg=get<0>(h);
									FIPSTER_ASSERT(hg->D ==1);
									if(grid1D) FIPSTER_ASSERT("not implemented");
									auto& val = f->additional.at(state_it.ind).second;
									res = val ? boost::lexical_cast<string>(hg->coords[0][val->ind]) 
										: "NaN"; 
								}else{
									auto& i = interpolated[j];
									res = !i ? "NaN" :
										boost::lexical_cast<string>(i->at(state_it.ind));
								}
								out<<res;
							}else
								out<<"--";
							double exercise = f->additional.at(state_it.ind).first ? 1 : 0;
							out<<"\t"<<exercise;
						}else
							out<<0;
						out<<"\t";
					}else{ //output hedging bands
						FIPSTER_ASSERT(bands2[j]);
						auto hg = get<0>(*(field[j].hedging));
						auto &bands = bands2[j]->at(state_it.ind);
						uint i=0,imax=2;
						for(; i<imax; i++)
							out<< boost::lexical_cast<string>(band_results[j][i].at(state_it.ind)) <<"\t";
						for(; i<4; i++)
							out<< ( bands[i] ? boost::lexical_cast<string>(hg->coords[0][bands[i]->ind])
											: "NaN" ) <<"\t";
						j+=hg->N[0];
					}
				}
				out<<"\n";
			}
			
			/*
			if(sol.output_as_geometric_average_grid_name){
				FIPSTER_ASSERT(grid.D==1);
				ofstream out2 (global_config::get().output_path+sol.field_name+"-geom_avg-"+
									*sol.output_as_geometric_average_grid_name+"-"+toS(time)+".csv");
				out2<<"#"<<sol<<endl;
				auto& grid2=*betterAt(sol.objs->grids,*sol.output_as_geometric_average_grid_name);
				grid_iterator<J_SGS,true> state_it(grid2);
				uint D=grid2.D;
				for(;state_it;++state_it){
					for(uint i=0;i<D;i++)
						out2<<boost::lexical_cast<string>(state_it.state()[i])<<"\t";
					out2<<boost::lexical_cast<string>(
						interpolate(grid.coords[0],**result,geometric_average(state_it.state(),D))
						)<<"\n";
				}
			}

			//result->print(out,grid);

			grid_iterator<J_SGS,true> state_it(sol.grid);
			for(;state_it;++state_it){
				//if(state_it.on_leading_boundary()) out<<"\n";
				for(uint i=0;i<grid.D;i++)
					out<<boost::lexical_cast<string>(state_it.state()[i])<<"\t";
				out<<boost::lexical_cast<string>(result->at(state_it.ind))<<"\n";
			}
		*/
			
			out.close();
			if(0)thread_logger()<<"OUTPUT | END"<<endl;
			return result_sptr_t(new int());
		}
	};
		
		// ostream& (ostream& out,const boost::optional<double>& v){
		// 	out << v ? boost::lexical_cast<string>(v) : "NaN";
		// 	return out;
		// }


	solution::solution(const ptree& pt,shared_const_objs_t objs)
		: id(betterGet<string>(pt,"<xmlattr>.id","solution"))
		,grid_name(betterGet<string>(pt,"grid",id))
		,fields_tree(betterChild(pt,"fields",id))
		,fixed_plot_time(getBtime(pt,"fixedPlotTime",id))
		,times(betterChild(pt,"times",id),id)
		,active(pt.get("<xmlattr>.active",true))
		,output_as_geometric_average_grid_name
		(pt.get_optional<string>("outAsGeometricAverageOn"))
		,show_hedging_pos(pt.get("<xmlattr>.hedgingPos",false))
		,objs(objs)
	{
		auto pt_single_values = pt.get_child_optional("singleValues");
		if(pt_single_values)
			single_values=read_vecs<double>("at","x",*pt_single_values
																			,"singleValues");
	}

	void solution::setup_nodes()
	{
		if(!active) return;

		// get grid (once)
		if(!grid) grid = betterAt(objs->grids,grid_name,id);
		
		//output gnuplot starting batch file:
		ofstream out (global_config::get().output_path+"plot"+
									toS(global_config::get().plot_count++)+".sh",ios_base::trunc);
			out<<"theta="<<times.start<<";stepsize="<<times.step_size<<
				";filename='"<<id<<"-"<<grid_name<<"';";
			if(fixed_plot_time)
				out<<"time='"<<*fixed_plot_time<<"'";
		out.close();

		read_referenced_fields(id,grid,fields_tree,fields,objs);

		//create vector of times
		not_done = times.times();
		
		//add fixed_plot_time if not already in the list of times
		if(fixed_plot_time)
			not_done.insert(*fixed_plot_time);

		//create nodes
		vector<spacetime_field<J_SGS>::sender_ptr> senders;
		for(auto it : not_done)
		{
			senders.clear();
			for(auto& j : fields)
				senders.push_back(j.get_node(it,grid));
			create_node(result_dumper_body(this,it),
									new_runtime_join_node(senders).get());
		}
			

		task_queue::finish_tasks();
	}
}

	ostream& operator<<(ostream& out,const solution& sol){
		//return out<<sol.name<<" at time "<<sol.time<<" on grid "<<sol.grid_name;
		for(uint i=1;i<=sol.grid->D;i++)
			out<<"sv"<<i<<"\t";
		for(auto& i : sol.fields){
			if(i.id){
				if(sol.show_hedging_pos && i.hedging && !get<2>(*i.hedging)){
					double h = get<0>(*i.hedging)->coords[0][get<1>(*i.hedging).ind];
					for(uint j=0; j<3; j++)
						out<<h<<"\t";
				}else
					out<<i.id<<"\t";
			}
		}
		return out;
	}

}


#include "type_traits.h"

