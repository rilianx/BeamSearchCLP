#include "SortedList.h"

 using namespace odp;
 namespace odp{
	 
	 void Standard_SortedList::update_rank(int nb_blocks, double w, map<BoxType*,int>& nb_left_boxes, Container &cont, FSpace* fsp_null) const{
   	    list<Block*>::const_iterator it=sorted_blocks.begin();
	    int i=0;
            for(;it!=sorted_blocks.end();it++){
	       if((*it)->ranking<0.0) {continue;}
	       (*it)->ranking+= w*(double(i)/double(nb_blocks));
	       i++;
	    }
	 }
	 
       void Vloss_SortedList::update_rank(int nb_blocks, double w, map<BoxType*,int>& nb_left_boxes, Container &cont, FSpace* fspace ) const{
	   if(knapsack) cont.solveKnapsack(nb_left_boxes);
    	   list<Block*>::const_iterator it=sorted_blocks.begin();
    	   long64 min_fitness=-1;
    	   long64 max_fitness= 0;
    	   for(;it!=sorted_blocks.end();it++){
    		  if((*it)->ranking<0.0){  continue; } //the block does not fit in the space
		  
    		 (*it)->fitness=cont.compute_VLossFitness(*it, *fspace);
    		 if((*it)->fitness>max_fitness) max_fitness=(*it)->fitness;
    		 if((*it)->fitness<min_fitness)
    			 min_fitness=(*it)->fitness;
		 

    	   }
	   
    	   long64 range_fitness=max_fitness-min_fitness;
    	   if(range_fitness==0) return;
	   
    	   for(it=sorted_blocks.begin();it!=sorted_blocks.end();it++){
    		if((*it)->ranking<0.0) continue;  
    		(*it)->ranking+= w*(double(max_fitness-(*it)->fitness)/double(range_fitness));
    		
    	   }
       }

       void actualize_rank(int nb_blocks, list<Block*>& blocks, vector<SortedList*>& lists, double* weights, 
	 map<BoxType*,int>& nb_left_boxes, Container& cont, FSpace* fspace){
	      list<Block*>::iterator it=blocks.begin();
	    	for(;it!=blocks.end();it++){
	    	      if((*it)->ranking<0.0) continue;	
                  else 	(*it)->ranking=0.0;	
		}
  	     for(int i=0; i<lists.size(); i++)
  		 if(weights[i]>0) lists[i]->update_rank(nb_blocks, weights[i], nb_left_boxes, cont, fspace);
  
             it=blocks.begin();
             //for(it=blocks.begin();it!=blocks.end();it++)
               // cout << (*it)->ranking << endl;
             //exit(0);
           
	 }	 
	 
	 
 }
