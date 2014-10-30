#ifndef _CONTAINER
#define _CONTAINER

#include "btBulletDynamicsCommon.h"
#include <BulletCollision/BroadphaseCollision/btDbvt.h>
#include "Policy.h"
#include "Box.h"
#include "Block.h"
#include <iostream>
#include <map>
#include <set>
#include <list>
#include <stack>

using namespace std;
namespace odp {




  
  struct compare{   
    bool operator() (const Box & sp1, const Box & sp2) const {
		if(sp1.getY()!=sp2.getY()) return (sp1.getY()<sp2.getY());
		if(sp1.getX()!=sp2.getX()) return (sp1.getX()<sp2.getX());
		if(sp1.getZ()!=sp2.getZ()) return (sp1.getZ()<sp2.getZ());

		if(sp1.Lengths().y()!=sp2.Lengths().y()) return (sp1.Lengths().y()<sp2.Lengths().y());
		if(sp1.Lengths().x()!=sp2.Lengths().x()) return (sp1.Lengths().x()<sp2.Lengths().x());
		return (sp1.Lengths().z()<sp2.Lengths().z());
			 
    }
  };
 
  struct by_manhattan_distance{  
    bool operator() (const FSpace & sp1, const FSpace & sp2) const{
       if(sp1==sp2) return false;
       if(sp1.manhattan_distance != sp2.manhattan_distance) return sp1.manhattan_distance < sp2.manhattan_distance;
			 
			 /* list<long64>::const_iterator it1=sp1.anchor_dist.begin();
       list<long64>::const_iterator it2=sp2.anchor_dist.begin();
	     if(*it1!=*it2) return (*it1<*it2);
	     it1++; it2++;
	     if(*it1!=*it2) return (*it1<*it2);
	     it1++; it2++;
	     if(*it1!=*it2) return (*it1<*it2);*/
	   	   
	    if(sp1.get_volume() != sp2.get_volume()) return (sp1.get_volume() > sp2.get_volume());
	   
       return compare()(sp1,sp2);

    }
 };  
 
  struct by_anchor_distance{  
    bool operator() (const FSpace & sp1, const FSpace & sp2) const{
       if(sp1==sp2) return false;
       
       list<long64>::const_iterator it1=sp1.anchor_dist.begin();
       list<long64>::const_iterator it2=sp2.anchor_dist.begin();
	   if(*it1!=*it2) return (*it1<*it2);
	   it1++; it2++;
	   if(*it1!=*it2) return (*it1<*it2);
	   it1++; it2++;
	   if(*it1!=*it2) return (*it1<*it2);
	   	   
	   if(sp1.get_volume() != sp2.get_volume()) return (sp1.get_volume() > sp2.get_volume());
	   
       return compare()(sp1,sp2);

    }
 };




  class Container{
    
  public:
    enum code{OK, COLLISION, ALREADY_IN, NO_GAP};
    
    static double alpha;
    static bool fsb;
    
    /********** gaps/boxes structure **************/
    btDbvt _tree_fspaces;
    
    set<FSpace, by_manhattan_distance> free_spaces;
    map<btDbvtNode *, set<FSpace, by_manhattan_distance>::iterator> node2space;


    btDbvt _tree_blocks;
    map<btDbvtNode*, Block_obj*> node2block;
    set<Block_obj*> removable_blocks;
    list<Block_obj*> locations;
        
    long64 L, W, H;
    
    //for the knapsack solutions
    int *mL, *mW, *mH;
    /****************************************/
    stack< list<FSpace> > rem_spaces;
    long64 _occupied_volume0;
    
    
    Container(int L, int W, int H);
      
    
	//Creates a new container.
	//All the blocks and free spaces are copied
    Container(Container& cont);


    //creates the object block and its supporting blocks
    Block_obj* create_bloxR(Block_obj* b, map<Block_obj*,Block_obj*> &created);
    
    void print_fspaces();
    void print_blocks();
    
    void remove_last_block(map<BoxType*,int>& nb_left_boxes); 


    //crea un nuevo objeto Block_obj (copia de blox) que sera asociado al contendor
    //remspaces=true indica que los espacios libres que no pueden contener ningun bloque pueden ser eliminados del contenedor
    int insert(Block_obj* blox, map<BoxType*,int>& nb_left_boxes, bool remspaces=true);  
    
	//retorna el volumen del espacio mas grande
    long64 remove_unfeasible_spaces(map<BoxType*,int>& nb_left_boxes);
    bool is_feasible(const Box& fsp, map<BoxType*,int>& nb_left_boxes);


    void remove(Block_obj* blox, map<BoxType*,int>& nb_left_boxes);
    void insert_spaces(list<Box>& initial, list<Box>& news,list<Box>::iterator &ref);
    void insert_spaces(list<Box>& initial, const Box& new_box,list<Box>::iterator &ref);

    void remove_fspace(const FSpace& fspace) {
         _tree_fspaces.remove(fspace.node);
         free_spaces.erase(fspace);
    }
    
    long64 occupied_volume(){
       long64 volume=0;

       list<Block_obj*>:: iterator it=locations.begin();
       for(;it!=locations.end();it++)
           volume+=(*it)->block->occupied_volume;
       
       return volume+_occupied_volume0;
    }
	 

    void solveKnapsack(map<BoxType*,int>& nb_left_boxes);
    long64 compute_VLossFitness(Block* block, const FSpace& fspace);
    long64 compute_VLossCorrected(Block* block, const FSpace& free_space);
    long64 compute_VLoss(Box& box, Box& free_space);

    void insert_fspace(FSpace &fspace){
          btDbvtNode * node = _tree_fspaces.insert(fspace, NULL);
  	    node2space[node]=free_spaces.insert(FSpace(fspace, node)).first;
    }
    

  
   ~Container(){
	 
	 while(!locations.empty()){
		 delete locations.front();
		 locations.pop_front();
	 }

   }
   
     static void remove_unfeasible_blocks(list<Block*>& blocks, map<BoxType*,int>& nb_left_boxes, set<FSpace, by_manhattan_distance>& fspaces);
     static void remove_unfeasible_blocks_conserv(list<Block*>& blocks, map<BoxType*,int>& nb_left_boxes, set<FSpace, by_manhattan_distance>& fspaces);


   bool fit_spaces(Block_obj& b1, Block_obj& b2);
   
  private:
    
    void block1(btDbvtAabbMm& _outbox, btDbvtAabbMm& _inbox, Box& blox, list<Box>& new_gaps, Block* block);
	
    
    void _remove_fspace(btDbvtNode* node); 
    

    
  };

  /** Marca bloques para eliminar e inicializa el resto **/
  void initialize_rankings(list<Block*>& blocks, map<BoxType*,int>& nb_left_boxes, set<FSpace, by_manhattan_distance>& fspaces);

  /** Elmina bloques marcados para eliminación de la lista**/
  void remove_unfeasible_blocks_rank(list<Block*>& blocks);

  void remove_unfeasible_blocks(list<Block*>& blocks, map<BoxType*,int>& nb_left_boxes);
  bool by_volume2(Box & b1, Box & b2);
}

#endif
