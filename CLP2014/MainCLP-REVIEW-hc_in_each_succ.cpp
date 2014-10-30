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
#include "SortedList.h"
#include "CLP.h"
#include <BulletCollision/BroadphaseCollision/btDbvt.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <time.h>
#include <vector>
#include <tr1/unordered_set>


using namespace std;
using namespace std::tr1;
using namespace odp;


void add(map<BoxType*,int> &boxes1, const map<BoxType*,int> &boxes2) {
	map<BoxType*,int>::const_iterator it_nb;
    for(it_nb = boxes2.begin(); it_nb!=boxes2.end(); it_nb++)
		boxes1[(*it_nb).first]+=(*it_nb).second;
};

bool sub(map<BoxType*,int> &boxes1, const map<BoxType*,int> &boxes2) {
	map<BoxType*,int>::const_iterator it_nb;
    for(it_nb = boxes2.begin(); it_nb!=boxes2.end(); it_nb++){
		boxes1[(*it_nb).first]-=(*it_nb).second;
		if(boxes1[(*it_nb).first]<0) return false;
		if(boxes1[(*it_nb).first]==0) boxes1.erase((*it_nb).first);
	}
	return true;
};

long64 compute_fitness(CLP& clp, Container& cont, map<BoxType*,int>& nb_left_boxes, list<Block*>& blocks, long64& BEST_VOLUME);
long64 greedy(
  list<Block*>& blocks, double *w, Container& cont, map<BoxType*,int>& nb_left_boxes);
void beam_search(CLP& clp, set<Block*, compareBlocks>& set_blocks,  int w, long64 &best_volume);
void fair_beam_search(CLP& clp, set<Block*, compareBlocks>& set_blocks,  int w, long64 &best_volume);


//Nodo del arbol de busqueda
class Node{
  public:
		
		static long64 L;
		static long64 W;
		static long64 H;

        Container* cont;
        map<BoxType*,int>* nb_left_boxes;
  		
	    Block_obj* blox;
	    Node* father;


        long64 fitness; //volume estimado (ub) 
        long64 vol; //volumen ocupado

		
	    list<Block*> blocks;
		
		
		Node(long64 fitness, long64 vol, Container* cont, map<BoxType*,int>* nb_left_boxes, Node* father) : 
        fitness(fitness), 
		vol(vol), cont(cont), nb_left_boxes(nb_left_boxes), father(father), blox(NULL) {	}
	
		
		Node(Node* father, Block_obj* blox, long64 vol, long64 fitness) : fitness(fitness), vol(vol), cont(NULL), nb_left_boxes(NULL), 
		      father(father), blox(blox){
		}
		
       ~Node(){
		  if(cont) {
			delete cont; 
			delete nb_left_boxes;
			//if(blox) delete blox;
		  }	
	}
		
		//root node
		Node(Container* cont, map<BoxType*,int>* nb_left_boxes, list<Block*> blocks) : cont(cont), 
		nb_left_boxes(nb_left_boxes), fitness(0), vol(0), father(NULL), blocks(blocks), blox(NULL){
			
		}
		
		void hard_copy(){
			cont=new Container(*father->cont);
			nb_left_boxes=new map<BoxType*,int>(*father->nb_left_boxes);
			blocks=father->blocks;		
		}
		
		void hard_uncopy(){
			cont=NULL;
			nb_left_boxes=NULL;
			blocks.clear();		
		}
		
		void solveKnapsack(){
			 cont->solveKnapsack(*(nb_left_boxes)); 
		}
		
		void insert_Block(Block_obj* blox){
		  vol+=blox->block->occupied_volume;
		  cont->insert(blox, *(nb_left_boxes), false);
	    }
	    
	    void remove_Block(Block_obj* blox){
			vol-=blox->block->occupied_volume;
			cont->remove(blox,*(nb_left_boxes));
		}

        //imposible volver al estado anterior usando remove_Block -> para usar con greedy
		void insertOPT_Block(Block_obj* blox){ 
		  cont->insert(blox, *(nb_left_boxes));
	    }
	    
	    long64 evaluate(CLP& clp, long64& best_volume){
			return compute_fitness(clp, *(cont), *(nb_left_boxes), blocks, best_volume);
		}
		
		void remove_unfeasible_blocks(){
			Container::remove_unfeasible_blocks_conserv(blocks, *nb_left_boxes, cont->free_spaces);
		}
		
		void swap(Block_obj& b1, Block_obj& b2){
			cont->remove(&b1,*(nb_left_boxes));
			cont->remove(&b2,*(nb_left_boxes));
	        Block_obj tmp(b1);
			b1=Block_obj(b1.block,&b2); 
			b2=Block_obj(b2.block,&tmp);
			insert_Block(&b1);	
			insert_Block(&b2);	
		}


};

	   long64 Node::W=0;
	   long64 Node::L=0;
	   long64 Node::H=0;




int greedy_calls=0;
double GR_PERC=0.9;
int gruns=0;
int beams=-1;
bool time25=false, time50=false, time100=false;
void get_bloxs(Node& node, list<Block_obj*>& bloxs, int w);

int main(int argc, char** argv)
{

	ifstream in(argv[1]);
	int _inst0=atoi(argv[2]); //instancia de partida (e.g. 25)
	int _inst1=atoi(argv[3]); //instancia de partida (e.g. 25)
	double min_fr=atof(argv[4]);
    Container::fsb = (atoi(argv[5])==1);
    beams=-1; //search effort
	int time=atoi(argv[6]);
	string output=argv[7];

 	long64 w0=10;
    srand(1);
     
	Container::alpha=1.0;
	
	string line;
	getline(in, line ); //number of instances



     for(int inst=0;inst<100; inst++){
	   CLP clp(in,time); //se inicializa el problema, 600 segundos de tiempo
	   Node::W=clp.W;
       Node::L=clp.L;
       Node::H=clp.H;
	   
		 double time0=clp.get_time();
		 time25=false;
		 time50=false;
		 time100=false;

	   if(_inst0!=-1 && (inst<_inst0 || inst >_inst1)) continue;	
	   long64 best_volume=0;
 
     long64 w=w0;
     set<Block*, compareBlocks> set_blocks = clp.general_block_generator(min_fr,10000,Container::fsb); 
     set<Block*, compareBlocks>::iterator it=set_blocks.begin();
	 for(int i=0;it!=set_blocks.end();i++,it++) (*it)->id=i;    	
 	


        cout << "instance " << inst << endl;
        cout << "curr_w	current_util	current_time" << endl;
		 while(!clp.timeout() && w <= 100000){
		    beam_search(clp, set_blocks, w, best_volume);
			cout << w << ":	"<< double(best_volume) / double(clp.W*clp.L*clp.H) << "	" << (clp.get_time()-time0) << endl;
		    if(!clp.timeout()){
				if(beams!=-1) w=w*2.0;
				else w=ceil(double(w)*sqrt(2));
				

			 } 

	   }
		    	      
     //cout << double(best_volume) << endl;
		 cout <<"Utilization:" << double(best_volume) / double(clp.W*clp.L*clp.H) << " " << (clp.get_time()-time0) << endl;
         
         
         while(set_blocks.size()>0){
   		     delete *set_blocks.begin();
   		     set_blocks.erase(set_blocks.begin());
         }
	   
      }
}

bool node_by_fitness(Node*  n1, Node*  n2){
		if(n1->fitness !=	n2->fitness ) return n1->fitness >	n2->fitness ;
		return (n1->vol < n2->vol);
    
  }

struct byfitness {
  bool operator() (const Node*  n1, const Node*  n2) const
  {
		if(n1->fitness !=	n2->fitness ) return n1->fitness >	n2->fitness ;
		return (n1->vol < n2->vol);
  }
};

//~ void fair_beam_search(CLP& clp, set<Block*, compareBlocks>& set_blocks,  int w, long64 &best_volume){
			//~ 
	    //~ Container* cont= new Container(clp.L,clp.W,clp.H);
	    //~ map<BoxType*,int>* nb_left_boxes= new map<BoxType*,int>(clp.nb_boxes); //cajas que faltan por colocar
	    //~ list<Block*> blocks(set_blocks.begin(), set_blocks.end());   //bloques
	    //~ 
		//~ long64 alpha=double(clp.L*clp.W*clp.H)/M;
		//~ 
		//~ set<Node*,byfitness> S[M]; 
		//~ 
	//~ 
		//~ //nodo raiz
		//~ Node* node=new Node(cont, nb_left_boxes, blocks);
			//~ 
		//~ S[0].push_back(node);
					//~ 
//~ 
        //~ while(true && !clp.timeout()){
		    //~ set<Node*,byfitness> SS;
			//~ for(int i=0;i<M;i++){
			  //~ for(set<Node*,byfitness>::iteratior it=S[i].begin();it!=S[i].end();it++){
				  //~ if(!visited[*it]) {SS=S[i];  goto end_loop;}
			  //~ }
		    //~ }
		    //~ break; //no new nodes in the sets
		    //~ end_loop:
			//~ 
			//~ 
			//~ for(set<Node*,byfitness>::iteratior it=SS.begin();it!=SS.end();it++){
				//~ if(visited[*it]) continue;
				//~ Node* node=*it;
				//~ visited[node]=true;
				//~ 
				//~ node->solveKnapsack(); 
				//~ list<Block_obj> bloxs;
				//~ 
				//~ get_bloxs(*node, bloxs, w); //the best block-space pairs are obtained and put in succ
				//~ 
				//~ for(;bs!=bloxs.end() && !clp.timeout();bs++){
					//~ node->insert_Block(*bs);
					//~ long64 vol=node->vol;
					//~ list<Block_obj> path;
					//~ long64 eval=node->evaluate(clp, best_volume, path); //debe retornar lista de bloques colocados
					//~ int lastj=-1;
					//~ Node* lastn=node;
					//~ for(int j=0;j<path.size()+1;i++){
						//~ if(!created[make_pair(vol,eval)] && eval > worst_eval(S[vol/alpha])){
							//~ Node* n=new Node(lastn); //whole-copy
							//~ for(int k=max(0,lastj); k<j; k++)
							   //~ n->insert_Block(path[k]);
							//~ n->fitness=eval;
							//~ n->remove_unfeasible_blocks();
							//~ 
							//~ S[vol/alpha].insert(n);
							//~ if(S[vol/alpha].size()>w) {
								//~ set<Node*>::iteratior itt=SS.end(); --itt; 
								//~ if(!visited[*itt]){
								  //~ visited[*itt]=true;
								  //~ (*itt)->delete_structures();
							    //~ }
								//~ removeNodes.push_back(*itt);
								//~ S[vol/alpha].erase(*itt);
							//~ }
							//~ 
							//~ lastn=n;
						//~ }
					//~ }
					//~ node->remove_lastBlock();
					//~ 
				//~ }
				//~ 
				//~ 
				//~ node->delete_structures();
			//~ }
			//~ 
		    //~ while(!removeNodes.empty()){
				//~ delete removeNodes.front();
				//~ removeNodes.pop_front();
		   //~ }	
			//~ 
		//~ }
//~ 
		//~ for(int i=0;i<M;i++){
			//~ for(set<Node*,byfitness>::iteratior it=S[i].begin();it!=S[i].end();it++)
				//~ delete *it;  
		//~ }		
			//~ 
//~ }

struct pair_hash {
    inline std::size_t operator()(const std::pair<int,int> & v) const {
        return v.first*31+v.second;
    }
    inline std::size_t operator()(const std::pair<ulong64,ulong64> & v) const {
        return v.first*31+v.second;
    }
};

bool hillclimbing(CLP& clp, Node* n, long64 &best_volume){
	vector<Block_obj*> removable_blocks(n->cont->removable_blocks.begin(),n->cont->removable_blocks.end());
	bool success=false;
	int m=removable_blocks.size();
	unordered_set< pair<int,int>, pair_hash > tabulist;

	if (m<=1) return false ;
	
	int tries = log((double) m)/log(2.0);
	
	for(int i=0; i<tries; i++){
		int b1=rand()%(m-1);
		int b2=m-1;//b1+1+rand()%(m-b1-1);
		bool flag=false;
		int b10=b1, b20=b2;
		while(!flag){
			// cout << b1 << "," << b2 << endl;
		  if((removable_blocks[b1]->block->w!=removable_blocks[b2]->block->w || 
		    removable_blocks[b1]->block->l!=removable_blocks[b2]->block->l || 
		    removable_blocks[b1]->block->h!=removable_blocks[b2]->block->h ) &&
		    tabulist.find(make_pair(b1,b2))==tabulist.end() &&
		    n->cont->fit_spaces(*removable_blocks[b1],*removable_blocks[b2])){
		       flag=true;  
		  }else{
		    b2++;
		    while(b2>=m){
		      b1=(b1+1)%(m-1);
		      b2=b1+1;
		    }
		    if(b10==b1 && b20==b2) return success; //no more valid movements
		  }
	    }

	    n->swap(*removable_blocks[b1],*removable_blocks[b2]);

	    ulong64 eval=n->evaluate(clp, best_volume);
	    
	    if(eval >= n->fitness){
			if(eval > n->fitness) i=0;
			n->fitness=eval;
			tabulist.clear();
			tabulist.insert(make_pair(b1,b2));
			success=true;

		}else{
			tabulist.insert(make_pair(b1,b2));
			n->swap(*removable_blocks[b1],*removable_blocks[b2]);

		}
			    			
	}
	return success;
}

void beam_search(CLP& clp, set<Block*, compareBlocks>& set_blocks,  int w, long64 &best_volume){
			
	    Container* cont= new Container(clp.L,clp.W,clp.H);
	    map<BoxType*,int>* nb_left_boxes= new map<BoxType*,int>(clp.nb_boxes); //cajas que faltan por colocar
	    list<Block*> blocks(set_blocks.begin(), set_blocks.end());   //bloques
	    unordered_set< pair<ulong64,ulong64>, pair_hash > volfit_nodes;
		
		list<Node*> S; 
	
		//nodo raiz
		Node* node=new Node(cont, nb_left_boxes, blocks);
			
		S.push_back(node);
					

		//S es el conjunto de nodos mantenidos en memoria
		while(S.size()!=0 && !clp.timeout()){

			list<Node*>::iterator itS=S.begin();		
			//N es usado para guardar los nodos de la siguiente generación
			list<Node*> N;
			//cout <<"spaces:" << (*itS)->cont->free_spaces.size() << endl;	
			//se expanden los nodos de la lista S
			for(;itS!=S.end();itS++){

				Node* node=*itS;
				list<Block_obj*> bloxs;
				
				node->solveKnapsack(); 
                

				get_bloxs(*node, bloxs, w); //the best block-space pairs are obtained and put in succ
         
         		list<Block_obj*>::iterator bs=bloxs.begin();
         		
         		Node* n_tmp=NULL; bool new_copy=true;
         		int i=0;
				for(;bs!=bloxs.end() && !clp.timeout();bs++){
                    if(new_copy){
						n_tmp=new Node(node, NULL, node->vol, node->fitness);
         		        n_tmp->hard_copy();
         		        new_copy=false;
					}
          
					n_tmp->insert_Block(*bs);
					n_tmp->fitness=n_tmp->evaluate(clp, best_volume);
					if(hillclimbing(clp, n_tmp, best_volume)){
						N.push_back(n_tmp);
						new_copy=true;
					}else{
					    Node* n=new Node(node, *bs, n_tmp->vol, n_tmp->fitness);
					    N.push_back(n);

					    n_tmp->remove_Block(*bs);
					    n_tmp->fitness=node->fitness;
					}
				}
				if(new_copy==false && n_tmp!=NULL) n_tmp->hard_uncopy();

					if(!time25 && clp.get_time()>30){
						time25=true;
						cout <<"Utilization(30):" << double(best_volume) / double(clp.W*clp.L*clp.H) << endl;
						
					}else if(!time50 && clp.get_time()>150){
						time50=true;
						cout <<"Utilization(150):" << double(best_volume) / double(clp.W*clp.L*clp.H) << endl;
					}else if(!time100 && clp.get_time()>300){
						time100=true;
						cout <<"Utilization(300):" << double(best_volume) / double(clp.W*clp.L*clp.H) << endl;
					} 

				
			}


			//the best successors are maintained
			N.sort(node_by_fitness); 
			list<Node*>::iterator itN=N.begin();

			long64 old_fitness;
			for(int i=0; itN!=N.end(); i++, itN++){
				//se mantienen solo los w mejores nodos, el resto se elimina
				if(i>= ((beams==-1)? w:beams)){
					//se sobrepaso el tamano de la lista de beams (se eliminan los nuevos nodos)
					delete (*itN)->blox;
					delete *itN;
					itN=N.erase(itN); itN--;
				}else{				
					Node *n = *itN;
					if(i>0 && volfit_nodes.find(make_pair(n->vol,n->fitness))!=volfit_nodes.end() /*n->fitness==old_fitness*/){ 
						//se eliminan nodos idénticos (o casi)
						delete (*itN)->blox;
						delete *itN; 
						itN=N.erase(itN); itN--; i--; 
					}else{
						//se genera el sucesor seleccionado, se copia el contenedor, 
						//las cajas restantes (tipo y cantidad) y la lista de blocks factibles
						if(n->blox){
						   n->hard_copy();
 					       n->insert_Block(n->blox);
					    }
						n->remove_unfeasible_blocks();
 					    old_fitness=n->fitness;
 					    volfit_nodes.insert(make_pair(n->vol,n->fitness));
					}
				}
				
			}
				

			while(!S.empty()){
				delete S.front();
				S.pop_front();
			}	

			S=N;
				
		} //hasta aqui el beam search
		
		//Se limpia la memoria ¿completamente?
		while(!S.empty()){
				delete S.front();
				S.pop_front();
		}				

			

}






   //retorna numero de bloques que calzan en el espacio correspondiente
   int actualize_rankings_fspace(list<Block*>& blocks, const FSpace& fspace, map<BoxType*,int>& nb_left_boxes){
   	int nb_blocks=0;
      list<Block*>::iterator it=blocks.begin();
   	for(;it!=blocks.end();it++){
   	      if((*it)->ranking==-2e10) continue;	
            if(!(*it)->feasible(nb_left_boxes)){
                   (*it)->ranking=-2e10;
      		 continue;		    
   		}else if(!(**it<=fspace))
   		    (*it)->ranking=-1e10;
   		else{
   	          (*it)->ranking=0;
   		    nb_blocks++;
   	      }		
      }
      return nb_blocks; 	
   }

   void selectBlocks_rank(list<Block*>& blocks, int n, map<BoxType*,int>& nb_left_boxes, 
	 set<Block*, by_ranking>& blocks_tmp){
	      
	blocks_tmp.clear();   
	
   	list<Block*>:: iterator it=blocks.begin();
   	for(;it!=blocks.end();it++){

	   if(((*it)->ranking) <= -1e10) continue;

   	   if(blocks_tmp.size()<n || by_ranking()(*it,(*blocks_tmp.rbegin())) ){

     	      if(blocks_tmp.size()==n){
       	     set<Block*, by_ranking>::iterator ii=blocks_tmp.end(); --ii;
     	           blocks_tmp.erase(ii);
     	      }
   	      blocks_tmp.insert(*it);

   	   }
      }
	      
   }
	 
 

   Block* selectBlock_greedyrank(list<Block*>& blocks,  int nb_blocks,
	  FSpace& fspace,  map<BoxType*,int>& nb_left_boxes){
	   
       
      set<Block*, by_ranking> blocks_tmp;
	  selectBlocks_rank(blocks, 1, nb_left_boxes, blocks_tmp);

	  set<Block*, by_ranking>::iterator it=blocks_tmp.end();
	  it--;
	  return *(it);
	  	  
   }
   
   Block* selectBlock_rank(list<Block*>& blocks,
	 FSpace& fspace,  map<BoxType*,int>& nb_left_boxes){

      return selectBlock_greedyrank(blocks, 0, fspace,  nb_left_boxes);
   }
  
  long64 greedy(
	  list<Block*>& blocks, Container& cont, 
        map<BoxType*,int>& nb_left_boxes){
	      
         while(!cont.free_spaces.empty()){
	
       	   FSpace fspace = *(cont.free_spaces.begin()); 
					 //K3: min anchor corn		   
		  
     	        //se coloca en el fspace la primera caja que quepa   
					 int nb_blocks=actualize_rankings_fspace(blocks, fspace, nb_left_boxes);
					 if(nb_blocks!=0){
						 
						 list<Block*>::iterator itB=blocks.begin();
						 for(;itB!=blocks.end();itB++){
							  if((*itB)->ranking<=-1e10) continue;
						 	  (*itB)->ranking=cont.compute_VLossFitness(*itB, fspace);
								 
								//1e14 + V - Vloss
						 }
						 
						 Block* block=selectBlock_greedyrank(blocks, nb_blocks, fspace, nb_left_boxes);
						 cont.insert(new Block_obj(block,fspace), nb_left_boxes);
						 
				    }else
						cont.remove_fspace(fspace);
		   
		   //cout << cont.occupied_volume() << endl;
        }
        return cont.occupied_volume();
     }   




     //create nodes with its associated action
		 void get_bloxs(Node& node, list<Block_obj*>& bloxs, int w){
			 initialize_rankings(node.blocks, *node.nb_left_boxes, node.cont->free_spaces);           
			 if(!node.father) w*=w;
			 
			 set<FSpace, by_manhattan_distance>::iterator itsp=node.cont->free_spaces.begin();	 
			 set<FSpace, by_manhattan_distance>::iterator sp;
			 if(itsp==node.cont->free_spaces.end()) return;		 

			 
			 int nb_blocks=0;	 int pos_space=0;
			 while(itsp!=node.cont->free_spaces.end() && nb_blocks<w){
			 	
			 			 list<Block*>::iterator itB=node.blocks.begin();
			 			 
			 			 for(;itB!=node.blocks.end();itB++){
			 					  //bloque no puede ser construido
			 				  if(((*itB)->ranking )<= -2e10) {continue;} 
			 						//bloque con ranking ya asignado
			 						
			 						if(pos_space!=0 && ((*itB)->ranking) != -1e10) continue; 
			 						
			 						
			 				 	  (*itB)->ranking=node.cont->compute_VLossFitness(*itB, *itsp); 
			 						if(((*itB)->ranking) != -1e10){
										
 		 						       nb_blocks++;
 		 						       
 		 						       
			 						   //(*itB)->pos_space=pos_space;
			 						   //(*itB)->location=(*itB)->get_location(*itsp);

			 						}
			 				 }
			 				 if(nb_blocks>0){ sp=itsp; break;}  //solo se consideran bloques que caben en el primer espacio factible
			 				 pos_space++;			
			 				 itsp++; 
                             
			 }
			 
   

			 set<Block*, by_ranking> blocks_tmp;
	     	 selectBlocks_rank(node.blocks, w, *node.nb_left_boxes, blocks_tmp);
			 
			 bool included_best_path=true; 
			 set<Block*, by_id>:: iterator it=blocks_tmp.begin();

			 for(; it!=blocks_tmp.end(); it++){	 
				 Block* b= *it;
                 bloxs.push_back(new Block_obj(b, *sp));	

			 }
																		
		 }

     long64 compute_fitness(CLP& clp, Container& cont, map<BoxType*,int>& nb_left_boxes,
		 list<Block*>& blocks, long64& BEST_VOLUME){
			 
			 
		   long64 max_fitness=cont.occupied_volume();
		   if(clp.timeout()) return max_fitness;
				 
		   Container cont_tmp(cont);
		   map<BoxType*,int> boxes_tmp=nb_left_boxes;
		   list<Block*> blocks_tmp=blocks;
		   cont_tmp.remove_unfeasible_spaces(boxes_tmp);
		   initialize_rankings(blocks_tmp, boxes_tmp, cont_tmp.free_spaces);
				 
		   max_fitness=greedy(blocks_tmp, cont_tmp, boxes_tmp);					 			 	 	 
				 
		  
		  
		  if(max_fitness>BEST_VOLUME){
			  BEST_VOLUME=max_fitness;
			 list<Block_obj*>::iterator it1=cont_tmp.locations.begin();

			 int collisions=0, outers=0;

			 while(it1!=cont_tmp.locations.end()){
			     list<Block_obj*>::iterator it2=it1; it2++;
				 while(it2!=cont_tmp.locations.end()){
					 if((*it1)->collides(**it2)) collisions++;
					 it2++;
				 }
				 if((*it1)->getX() < 0 || (*it1)->getXmax() > clp.L || (*it1)->getY() < 0 || (*it1)->getYmax() > clp.W || 
				    (*it1)->getZ() < 0 || (*it1)->getZmax() > clp.H ) outers++;
				 it1++;
			 }
			 if(collisions+outers >0)
			   cout << "(fatal error) collisions:" << collisions << ", out blocks:" << outers << endl;
			 else cout << "volume utilization (%):" << (cont_tmp.occupied_volume()) / double(clp.W*clp.L*clp.H) 
			 << "(" << (cont.occupied_volume()) / double(clp.W*clp.L*clp.H) << ")" 
			 << ", nb blocks:" <<  cont_tmp.locations.size() << "(" << cont.locations.size() << ")" 
			 << ", rem blocks:" <<  cont_tmp.removable_blocks.size() << "(" << cont.removable_blocks.size() << ")" 
			 << endl;
		  }
			 
		  return max_fitness;
     }



	
  




