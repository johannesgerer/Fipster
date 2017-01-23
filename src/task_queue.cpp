#include "precompiled.h"
#include "task_queue.h"


using namespace fipster::_task_queue;

void task_queue::finish_tasks(){
	auto& q = get().q;
	while(!q.empty()){
		// cout<<"q.size(): "<<q.size()<<endl;
		q.front()();
		q.pop();
	}
}
