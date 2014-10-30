/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2007 Erwin Coumans  http://continuousphysics.com/Bullet/

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#include "btBulletDynamicsCommon.h"
#include "Policy.h"
#include "Container.h"
#include "Block.h"
#include "CLP.h"
#include <BulletCollision/BroadphaseCollision/btDbvt.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <time.h>



using namespace std;
using namespace odp;

class Node{
   public:
   Node(Block *b, long64 x, long64 y, long64 z): block(b), x(x), y(y), z(z) {}
   long64 x;
   long64 y;
   long64 z;
   Block* block;
};

long64 greedy(list<Block*>& blocks, Container& cont, map<BoxType*,int>& nb_left_boxes, 
	double grasp_perc=1.0, list<Node>* path=NULL);



int main(int argc, char** argv)
{


	ifstream in(argv[1]);
	int _inst=atoi(argv[2]);
	
	string line;
	getline(in, line ); //number of instances

      for(int inst=0;inst<100; inst++){
	   CLP clp(in,600); //se inicializa el problema, 600 segundos de tiempo
	   if(_inst!= -1 && inst!=_inst) continue;	
	      
	    //Generación simple de bloques
	   //set<Block*, compareBlocks> set_blocks=clp.single_box_block_generator();
	   set<Block*, compareBlocks> set_blocks = clp.general_block_generator();

	   list<Block*> blocks(set_blocks.begin(), set_blocks.end());
         blocks.sort(by_volume); //necesario para heuristica que elije bloques
	  	
	   Container cont(clp.L,clp.W,clp.H);
         map<BoxType*,int> nb_left_boxes=clp.nb_boxes;	

         greedy(blocks, cont, nb_left_boxes);
	

         cout <<  (double(cont.occupied_volume()) / double(clp.W*clp.L*clp.H)) << endl;
		  
    }

}

//Selecciona el mejor bloque usando evaluación:  V - Vloss. V es el volumen del bloque
//Vloss, es el volumen que se perdería en el espacio al colocarlo
 Block* selectBlock_Vloss(list<Block*>& blocks, Container& cont, FSpace& fspace,  map<BoxType*,int>& nb_left_boxes, 
   double grasp_perc=1.0){
      long64 max_fitness=0;
      	
      Block* block=NULL;
      list<Block*>:: iterator it=blocks.begin();
     	for(;it!=blocks.end();){
   	       if(max_fitness >= 1e14 + (*it)->occupied_volume ) break;
           	
 	       if(!(*it)->feasible(nb_left_boxes)){
 	         it=blocks.erase(it);
 			 continue;
 		   }
		   
	 	   
	       long64 fitness=cont.compute_VLossFitness(*it, fspace);
     	      //long64 fitness=cont.compute_defrag(*it, fspace);
   	       if(fitness > max_fitness){ 
			   if(!block) block=*it; //bloque por defecto (max_volume)
			   if(grasp_perc < 1.0 && (double)rand()/(double)RAND_MAX > grasp_perc){it++; continue;}		   		   
			   max_fitness=fitness; 
			   block=*it;
		   }
		   it++;
   	    }
		return block;

   }
   


//Coloca bloque más grande en min anchor corner space
long64 greedy(list<Block*>& blocks, Container& cont, map<BoxType*,int>& nb_left_boxes, double grasp_perc, list<Node>* path){
    
    cont.solveKnapsack(nb_left_boxes);
    while(!cont.free_spaces.empty()){
  	   FSpace fspace = *(cont.free_spaces.begin()); //K3: min anchor corner
	   bool is_ok=false;

	   
	   //se coloca en el fspace la primera caja que quepa   
	   Block* block=selectBlock_Vloss(blocks, cont, fspace,  nb_left_boxes, grasp_perc);
	   
  	   if(block){
  	     cont.insert(block, fspace, nb_left_boxes, false);
		 if(path){
                path->push_back(Node(block, cont.locations.front().getX(),
                cont.locations.front().getY(), cont.locations.front().getZ()));
		 }
       }else
         cont.remove_fspace(fspace);


   }
   return cont.occupied_volume();
}







