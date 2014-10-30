 #ifndef _SORTLIST
 #define _SORTLIST

#include <iostream>
#include <list>
#include <vector>
#include "Block.h"
#include "Container.h"

using namespace std;
 namespace odp {
 
 
     class SortedList{
     public:
    	 list<Block*> sorted_blocks;
    	 bool (*compare)(Block* & b1, Block* & b2);
       
	 virtual void init(list<Block*>& blocks) = 0;
	 virtual void update_rank(int nb_blocks, double w, map<BoxType*,int>& nb_left_boxes, Container &cont, FSpace* fspace=NULL) const=0;    
	     
     };
 
 
     class Standard_SortedList : public SortedList{
       public:
       Standard_SortedList(bool (*compare)(Block* & b1, Block* & b2))   {
		 this->compare=compare;
	 }
	 
	 void init(list<Block*>& blocks){
		sorted_blocks=blocks; 
		sorted_blocks.sort(*compare);
	 }
	 
	 virtual void update_rank(int nb_blocks, double w, map<BoxType*,int>& nb_left_boxes, Container &cont, FSpace* fspace=NULL) const;
	 
	 
     };
     

     class Vloss_SortedList : public SortedList{
     public:
     	  Vloss_SortedList(bool knapsack=false) : knapsack(knapsack){
     	  	//sorted_blocks.sort(by_volume);
     	  } 

 	 void init(list<Block*>& blocks){
 		sorted_blocks=blocks; 
 	 }
	  
	  virtual void update_rank(int nb_blocks, double w, map<BoxType*,int>& nb_left_boxes, Container &cont, FSpace* fspace=NULL) const;
	  bool knapsack;
     };
     
     void actualize_rank(int nb_blocks, list<Block*>& blocks, vector<SortedList*>& lists, double* weights, 
	     map<BoxType*,int>& nb_left_boxes, Container& cont, FSpace* fspace);
     
 
 }
 #endif