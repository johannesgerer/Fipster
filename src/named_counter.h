#pragma once
#ifndef AUTOMATIC_PRECOMPILATION
#include <map>
#include <string>
#endif



namespace fipster { 

	using namespace std;

	struct named_counter{

		
		//Singleton provider
		static named_counter& get_singleton(){
			static named_counter singleton;
			return singleton;
		}

		map<string,int> counters;

		template<class T> static void increment(T&& s,int n=1){ 
			get_singleton().counters[forward<T>(s)]+=n;
		}

		template<class T> static void decrement(T&& s,int n=1){ 
			get_singleton().counters[forward<T>(s)]-=n;
		}

		template<class T> static int& get(T&& s){ 
			return get_singleton().counters[forward<T>(s)];
		}

		
		template<class T> static void reset(T&& s){ 
			get_singleton().counters[forward<T>(s)]=0;
		}


		 static void print();

	};

}
