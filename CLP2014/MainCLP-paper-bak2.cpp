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



using namespace std;
using namespace odp;

class Action{
   public:
		 Action(): x(0), y(0), z(0), block(NULL) {}  //for the root node
     Action(Block *b, long64 x, long64 y, long64 z): block(b), x(x), y(y), z(z) {}
     long64 x;
     long64 y;
     long64 z;
     Block* block;
};

class Node{
  public:
		
		Node(list<Action>& path, long64 fitness, long64 vol, Container* cont, map<BoxType*,int>* nb_left_boxes, Node* father) : 
		path(path), fitness(fitness), 
		vol(vol), cont(cont), nb_left_boxes(nb_left_boxes), father(father), best_vol(-1){
			
		}
		
		Node(Node* father, Action act, long64 vol) : fitness(vol), vol(vol), cont(NULL), nb_left_boxes(NULL), 
		father(father), action(act), best_vol(-1){
			
		}		
		
		//root node
		Node(Container* cont, map<BoxType*,int>* nb_left_boxes, list<Block*> blocks) : cont(cont), 
		nb_left_boxes(nb_left_boxes), fitness(0), vol(0), depth(0), father(NULL), blocks(blocks), best_vol(-1){
			
		}
		
		void instantiate(Container* cont2, map<BoxType*,int>* nb_left_boxes2){
			cont=cont2;
			nb_left_boxes=nb_left_boxes2;			
		}
		
		
		
    Container* cont;
    map<BoxType*,int>* nb_left_boxes;
		
		Action action;
		Node* father;
    list<Action> path; //desde el nodo hasta destino
    long64 fitness; //volume estimado (ub) 
    long64 vol; //volumen ocupado
		long64 best_vol; //mejor volumen que se alcanzó a través del nodo
		int depth;
		
		list<Block*> blocks;
};

class vertex{
public:
	
	double* x;
	int size;
	double eval;
	
	~vertex(){
		delete[] x;
	}
	
	vertex& operator=( const vertex& vv ) {

		if( this != &vv ) {
			delete [] x;
			size=vv.size; eval=vv.eval;
			x = new double[ size ];
			for( int i = 0; i < size; i++ )
				x[i] = vv[i];
		}

		return *this;
	}
	
	vertex(int i) : size(i), eval(-1.0){
		x=new double[i];
		for( int i = 0; i < size; i++ ) x[i]=0.0;
		
	}
	
	vertex(const vertex& v) : eval(v.eval), size(v.size){
		x=new double[size];
		for(int i=0; i<size; i++)
			x[i]=v[i];
	}
	
	double& operator[](int i) const{
		return x[i];
	}
	
};

void expand(Node& node, list<Node*>& successors, int w, vector<vertex>& weights, vector<SortedList*>& sorted_lists);

long64 compute_fitness(CLP& clp, Container& cont, map<BoxType*,int>& nb_left_boxes, list<Action>& path,
		 int w, int depth, int max_depth, vector<vertex>& weights, list<Block*>& blocks, vector<SortedList*>& sorted_lists);

long64 greedy(list<Block*>& blocks, Container& cont, map<BoxType*,int>& nb_left_boxes, 
	double grasp_perc=1.0, list<Action>* path=NULL);
long64 greedy(
  list<Block*>& blocks, vector<SortedList*>& sorted_lists, double *w, Container& cont, map<BoxType*,int>& nb_left_boxes, list<Action>* path, double grasp=0.0);
void beam_search(CLP& clp, set<Block*, compareBlocks>& set_blocks, vector<SortedList*>& sorted_lists, 
    vector<vertex>& weights, int w, long64 &best_volume, double grasp, list<Action>& best_path, 
			bool parallel_search);


int greedy_calls=0;
double GR_PERC=0.9;
int gruns=0;

bool G2LA=false;
//double beta1, beta0;
int ss; //numero de grupos
int main(int argc, char** argv)
{


	ifstream in(argv[1]);
	int _inst0=atoi(argv[2]); //instancia de partida (e.g. 25)
	int _inst1=atoi(argv[3]); //instancia de partida (e.g. 25)
	long64 w0=atoi(argv[4]);
  srand(atoi(argv[5]));
	//bool parallel=(atoi(argv[6])==1);
	ss=atoi(argv[6]);
	int time=atoi(argv[7]);


	Container::alpha=1.0;

	
	string line;
	getline(in, line ); //number of instances


	vector<SortedList*> sorted_lists;
	sorted_lists.push_back(new Standard_SortedList(by_volume));
	sorted_lists.push_back(new Standard_SortedList(by_sqr_volumes));
	sorted_lists.push_back(new Standard_SortedList(by_volume_div_nb_box));
	Vloss_SortedList* vlsl=new Vloss_SortedList(false);
	sorted_lists.push_back(vlsl);

	    
     for(int inst=0;inst<100; inst++){
	   CLP clp(in,time); //se inicializa el problema, 600 segundos de tiempo
		 double time0=clp.get_time();
	   if(_inst0!=-1 && (inst<_inst0 || inst >_inst1)) continue;	
	   long64 best_volume=0;
     list<Action> best_path;	   
     long64 w=w0;
     set<Block*, compareBlocks> set_blocks = clp.general_block_generator();
     set<Block*, compareBlocks>::iterator it=set_blocks.begin();
	   for(int i=0;it!=set_blocks.end();i++,it++) (*it)->id=i;    	
	   	
	   	  
	
	   vector<vertex> weights;
	    //weights.push_back( nelder_mead(clp, set_blocks, sorted_lists, 2, 0.01, best_path, best_volume) );
	   // nelder_mead(clp, set_blocks, sorted_lists, 2, 0.01, best_path, best_volume);
	   
	   vertex v0(4);
	   v0[0]=0; v0[1]=0; v0[2]=0; v0[3]=1.0;
	   weights.push_back(v0);
	 //  v0[0]=0; v0[1]=1; v0[2]=0; v0[3]=0;
	 //  weights.push_back(v0);
	 //  v0[0]=0; v0[1]=0.3; v0[2]=0.3; v0[3]=0.3;
	 //  weights.push_back(v0);
	   // v0[0]=0; v0[1]=0; v0[2]=0; v0[3]=1.0;
	   // weights.push_back(v0);	   

		 while(!clp.timeout() && w <= 100000){
		    beam_search(clp, set_blocks, sorted_lists, 
		     weights, w, best_volume, 0, best_path, 0);
					cout << w << ":"<< double(best_volume) / double(clp.W*clp.L*clp.H) << " " << (clp.get_time()-time0) << endl;
				if(!clp.timeout()) w=ceil(double(w)*sqrt(2));	

	   }
		    	      
     //cout << double(best_volume) << endl;
		 cout <<"Utilization:" << double(best_volume) / double(clp.W*clp.L*clp.H) << " " << (clp.get_time()-time0) << endl;
         while(set_blocks.size()>0){
   		     delete *set_blocks.begin();
   		     set_blocks.erase(set_blocks.begin());
         }
	   
      }




}




bool node_by_volume(Node*  n1, Node*  n2){
	if(n1->vol != n2->vol) return 	(n1->vol  >	n2->vol );
		return (n1->fitness > n2->fitness);
    
  }


bool node_by_fitness(Node*  n1, Node*  n2){
		if(n1->fitness != n2->fitness) return 	(n1->fitness  >	n2->fitness );
		return (n1->vol < n2->vol);
    
  }


void beam_search(CLP& clp, set<Block*, compareBlocks>& set_blocks, vector<SortedList*>& sorted_lists, 
    vector<vertex>& weights, int w, long64 &best_volume, double grasp, list<Action>& best_path, 
    bool parallel_search){
			
			int node_limit=w;
			
	    Container* cont= new Container(clp.L,clp.W,clp.H);
	    map<BoxType*,int>* nb_left_boxes= new map<BoxType*,int>(clp.nb_boxes);
	    list<Block*> blocks(set_blocks.begin(), set_blocks.end());
			for(int i=0; i<sorted_lists.size(); i++)  sorted_lists[i]->init(blocks);    
			
			list<Node*> S; 
			list<Node*> allNodes; //para segunda etapa->mejora
			
			Node* node=new Node(cont, nb_left_boxes, blocks);
      list<Action> new_path;
			node->cont->solveKnapsack(*(node->nb_left_boxes)); 	
			long64 first_fitness=compute_fitness(clp, *(node->cont), *(node->nb_left_boxes), 
				new_path, w, 0, (G2LA)? 1:0, weights, node->blocks, sorted_lists);
			//best_volume=first_fitness;
			S.push_back(node);
					
			int depth=1;
			while(S.size()!=0 && !clp.timeout()){

				list<Node*>::iterator itS=S.begin();
				
				list<Node*> N;
				for(;itS!=S.end();itS++){

					Node* node=*itS;
					list<Node*> succ;
				 

					node->cont->solveKnapsack(*(node->nb_left_boxes)); 		  
					expand(*node, succ, w, weights, sorted_lists);
          
					if(succ.size()==0){//terminal node!
						Node* father=node->father;
						while(father){
							if(father->best_vol < node->vol) father->best_vol=node->vol;
							father=father->father;
						}
					}

					list<Node*>::iterator it_s=succ.begin();
					for(;it_s!=succ.end();it_s++){
						Node* s=*it_s;
						node->cont->insert(s->action.block, s->action.x, s->action.y, 
							s->action.z, *(node->nb_left_boxes));	
							
						//calcula volumen alcanzable desde el nodo
	  
						list<Action> new_path;				
						long64 new_fitness=compute_fitness(clp, *(node->cont), *(node->nb_left_boxes), 
							new_path, w, 0, (G2LA)? 1:0, weights, node->blocks, sorted_lists);


						if(new_fitness > s->fitness){
							
							s->path=new_path;
							s->fitness=new_fitness;
							
						}
						
						if(new_fitness > best_volume) best_volume=new_fitness;
						
						
						node->cont->remove_last_block( *(node->nb_left_boxes));	
						
					}
					
					if(parallel_search && succ.size()>0){
					   succ.sort(node_by_fitness);
						 N.push_back(succ.front());
						 succ.pop_front();
						 while(!succ.empty()){
							 delete succ.front();
						 	 succ.pop_front();
						 }
					}else{
						 N.insert(N.begin(), succ.begin(), succ.end());
					 }
				}
				
				if(N.size()>0){
					N.sort(node_by_volume);

          //se eliminan nodos idénticos (o casi)
          list<Node*>::iterator itN=N.begin();
					long64 old_fitness=-1;
					long64 old_vol=-1;
			  	for(; itN!=N.end(); itN++){
					  Node *n = *itN;
					  if(n->vol==old_vol && n->fitness==old_fitness){ 	      
						  delete *itN; 
						  itN=N.erase(itN); itN--;
					  }
						old_fitness=n->fitness;
						old_vol=n->vol;
					}

				
					int s=1; //min(w/ss,int(N.size()));
					list<Node*>* sublist=new list<Node*>[s];
					double subsize=double(N.size())/double(s);
		
					list<Node*>::iterator fin=N.begin(); int i=0;
					for(int l=0; l<s; l++){
						list<Node*>::iterator ini=fin;
						for(;i<int(subsize*double(l+1)) && fin!=N.end(); fin++, i++);
						if(ini!=fin)
							sublist[l].insert(sublist[l].begin(), ini, (l==s-1)? N.end():fin );
						sublist[l].sort(node_by_fitness);
					}
				

					for(int l=0; l<s; l++){

						list<Node*>::iterator itN=sublist[l].begin();
						long64 old_vol, old_fitness; 
						int ww=ceil((double)w/(double)s);
				
						for(int i=0; itN!=sublist[l].end(); i++, itN++){
							if(i>= ((G2LA)?1:ww)){
								delete *itN;
								itN=sublist[l].erase(itN); itN--;
							}else{
											
								Node *n = *itN;

								Container* cont_tmp = new Container(*n->father->cont);
								map<BoxType*,int>* boxes_tmp=new map<BoxType*,int>(*n->father->nb_left_boxes);
							
						
								cont_tmp->insert(n->action.block, n->action.x, n->action.y, 
									n->action.z, *boxes_tmp, false);
								cont_tmp->remove_unfeasible_spaces(*boxes_tmp);
						
								n->instantiate(cont_tmp, boxes_tmp);
								n->depth=depth;
								n->blocks=n->father->blocks;
								remove_unfeasible_blocks(n->blocks, *boxes_tmp, cont_tmp->free_spaces);
						
								old_fitness=n->fitness;
								old_vol=n->vol;
							}

						}
					}
					
					N=sublist[0];
					for(int l=1; l<s ;l++){
						N.insert(N.end(), sublist[l].begin(), sublist[l].end());
					}
				}
			

				
				//no se puede eliminar ya que cont es copiado a los sucesores
				while(!S.empty()){
					if((S.front())->cont) {
						delete (S.front())->cont; delete (S.front())->nb_left_boxes;
					  (S.front())->cont = NULL; (S.front())->nb_left_boxes = NULL;
					}			
					allNodes.push_back(S.front());			
					//delete S.front();
					S.pop_front();
				}	

				S=N;
				
			}
			
			if(clp.timeout()){
				while(!S.empty()){
					if((S.front())->cont) {
						delete (S.front())->cont; delete (S.front())->nb_left_boxes;
					  (S.front())->cont = NULL; (S.front())->nb_left_boxes = NULL;
					}			
					allNodes.push_back(S.front());			
					//delete S.front();
					S.pop_front();
				}				
			}

/*
			long64 x0=0, y0=(best_volume-first_fitness)/5, x1=0, y1=0;	
			list<Node*>::iterator it; bool done=false;
			for (it=allNodes.begin(); it!=allNodes.end(); ++it){
				Node* n=*it;
				if(!done && n->best_vol!=-1 && (n->best_vol - n->fitness)==0 && n->fitness>0){	
					x1=n->vol;
					y1=(n->best_vol - n->fitness);
				  done=true;
				}
        delete n;					
			}
			
			if(done){
			  beta1=double(y1-y0)/double(x1-x0);
			  beta0=y1-beta1*x1;
				if(beta1>0.0) {beta1=0.0; beta0=0.0;}
			
			}
*/			  
			 
		  	
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

   void selectBlocks_rank(list<Block*>& blocks, vector<SortedList*>& sorted_lists, Container& cont, 
  int n, map<BoxType*,int>& nb_left_boxes, set<Block*, by_ranking>& blocks_tmp){
	      
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

   Block* selectBlock_greedyrank(list<Block*>& blocks, int nb_blocks, vector<SortedList*>& sorted_lists, 
   Container& cont, FSpace& fspace,  map<BoxType*,int>& nb_left_boxes, double perc){
	   
	  
	  int size=max(int(perc*double(nb_blocks)+0.9999),1);
	  int r=rand()%size;

        set<Block*, by_ranking> blocks_tmp;
	  selectBlocks_rank(blocks, sorted_lists, cont, 
	     r+1, nb_left_boxes, blocks_tmp);

	  set<Block*, by_ranking>::iterator it=blocks_tmp.end();
	  it--;
	  return *(it);
	  	  
   }
   
   Block* selectBlock_rank(list<Block*>& blocks, int nb_blocks, vector<SortedList*>& sorted_lists, 
   Container& cont, FSpace& fspace,  map<BoxType*,int>& nb_left_boxes){

	   return selectBlock_greedyrank(blocks, nb_blocks, sorted_lists, 
	      cont, fspace,  nb_left_boxes, 0);
   }
  
  long64 greedy(
	  list<Block*>& blocks, vector<SortedList*>& sorted_lists, Container& cont, 
        map<BoxType*,int>& nb_left_boxes, list<Action>* path, double grasp){
	      
         while(!cont.free_spaces.empty()){
	
       	   FSpace fspace = *(cont.free_spaces.begin()); //K3: min anchor corn		   
		  
     	        //se coloca en el fspace la primera caja que quepa   
					 int nb_blocks=actualize_rankings_fspace(blocks, fspace, nb_left_boxes);
					 if(nb_blocks!=0){
						 
						 list<Block*>::iterator itB=blocks.begin();
						 for(;itB!=blocks.end();itB++){
							  if((*itB)->ranking<=-1e10) continue;
						 	  (*itB)->ranking=cont.compute_VLossFitness(*itB, fspace);
								 
								//1e14 + V - Vloss
						 }
						 
						 //long64_rank(nb_blocks, blocks, sorted_lists, w,
						 // nb_left_boxes, cont, &fspace);
						 Block* block=selectBlock_greedyrank(blocks, nb_blocks, sorted_lists, cont, fspace, nb_left_boxes, grasp);

						 cont.insert(block, fspace, nb_left_boxes, false);
						 
						 if(path){
							 path->push_back(Action(block, cont.locations.front().getX(),
								 cont.locations.front().getY(), cont.locations.front().getZ()));
						 }
						 }else
							 cont.remove_fspace(fspace);
		   
		   //cout << cont.occupied_volume() << endl;
        }
        return cont.occupied_volume();
     }   




     //create nodes with its associated action
		 void expand(Node& node, list<Node*>& successors, int w, vector<vertex>& weights, vector<SortedList*>& sorted_lists){
				


			 initialize_rankings(node.blocks, *node.nb_left_boxes, node.cont->free_spaces);
			 //nb_blocks=actualize_rankings_fspace(blocks, *itsp, *node.nb_left_boxes);                


			 
			 set<FSpace, by_anchor_distance>::iterator itsp=node.cont->free_spaces.begin();	 
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
			 							 (*itB)->pos_space=pos_space;
			 							 (*itB)->location=(*itB)->get_location(*itsp);
			 						}
			 				 }
			 				 pos_space++;			
			 				 itsp++; 

							 break;  //solo se consideran bloques que caben en el 
							 //primer espacio factible
			 }
			 
   

			 set<Block*, by_ranking> blocks_tmp;
			 
			 selectBlocks_rank(node.blocks, sorted_lists, *node.cont, w, *node.nb_left_boxes, blocks_tmp);
       //cout << blocks_tmp.size() << endl;
			 bool included_best_path=false;
			 set<Block*, by_id>:: iterator it=blocks_tmp.begin();
			 for(; it!=blocks_tmp.end(); it++){
				 
				 Block* b= *it;
         //cout << b->pos_space << endl;

          //nodo: padre + action
				 Node* n;	
				 if(node.path.size()>0 && b==node.path.front().block){ 
					 //acción corresponde a la del best_path
					 n=new Node(&node, node.path.front(), node.vol+b->occupied_volume);
					 list<Action> path=node.path; 
					 path.pop_front();
					 n->path=path; n->fitness=node.fitness; //mantiene beth_path
					 included_best_path=true;
				 }else{
						//lugar donde será colocado el bloque (usando metodo 
					  //anchor_distance)
					  Box blox = b->location; 
						//Box blox = b->get_location(*itsp);
						
					  n= new Node(&node, Action(b, blox.getX(), blox.getY(), blox.getZ()), node.vol+b->occupied_volume);
					}	
				 
				 successors.push_back(n);	
			 }

			 if(!included_best_path && node.path.size()>0){
				 Node* n = new Node(&node, node.path.front(), node.vol+node.path.front().block->occupied_volume );
				 list<Action> path=node.path; 
				 path.pop_front();
				 n->path=path; n->fitness=node.fitness;	 //mantiene beth_path
				 successors.push_back(n);	
			 }
																		
		 }

     long64 compute_fitness(CLP& clp, Container& cont, map<BoxType*,int>& nb_left_boxes, list<Action>& path,
		 int w, int depth, int max_depth, vector<vertex>& weights, list<Block*>& blocks, vector<SortedList*>& sorted_lists){
			 
			 
		   long64 max_fitness=cont.occupied_volume();
			 if(clp.timeout()) return max_fitness;
				 
			 if(depth==max_depth){
				 Container cont_tmp(cont);
				 map<BoxType*,int> boxes_tmp=nb_left_boxes;
				 initialize_rankings(blocks, boxes_tmp, cont_tmp.free_spaces);
				 max_fitness=greedy(blocks, sorted_lists, cont_tmp, boxes_tmp, &path, 0);					 			 	 	 
			 }else{
		  	 set<FSpace, by_anchor_distance>::iterator itsp=cont.free_spaces.begin();	
		  	 int nb_blocks=0;     
		     while(nb_blocks==0 && itsp!=cont.free_spaces.end()){
		  		     initialize_rankings(blocks, nb_left_boxes, cont.free_spaces);
		  		     nb_blocks=actualize_rankings_fspace(blocks, *itsp, nb_left_boxes);
		  		     if(nb_blocks==0)
		 					    itsp++;	 	           
		  	 }
				 if(nb_blocks==0) return max_fitness;	 
				 
				 FSpace fspace = (*itsp);
		 		 set<Block*, by_ranking> blocks_tmp; //bloques en el nivel depth del árbol
		     //calcula ranking de las cajas restantes
				 list<Block*>::iterator itB=blocks.begin();
				 for(;itB!=blocks.end();itB++){
					  if((*itB)->ranking<=-1e10) continue;
				 	  (*itB)->ranking=cont.compute_VLossFitness(*itB, fspace);
				 }
		// actualize_rank(nb_blocks, blocks, sorted_lists, weights[0].x, nb_left_boxes, cont, &fspace);
		 		 selectBlocks_rank(blocks, sorted_lists, cont, w, nb_left_boxes, blocks_tmp);
		 		 		 
				 set<Block*, by_ranking>:: iterator it=blocks_tmp.begin();
				 for(; it!=blocks_tmp.end(); it++){						 
					 cont.insert(*it, fspace, nb_left_boxes);
					 
					 list<Action> new_path;
					 long64 fitness=compute_fitness(clp, cont, nb_left_boxes, new_path,
					 		 w, depth+1, max_depth, weights, blocks, sorted_lists);
					 if(fitness>max_fitness){
						 path=new_path;
		 				 path.push_front(Action(*it, cont.locations.front().getX(),	
							 cont.locations.front().getY(), cont.locations.front().getZ()));	 
						 max_fitness=fitness;
					 }
					 
					 cont.remove_last_block(nb_left_boxes);			 
				 }				 	
			 }
			 
			 return max_fitness;
     }



	
  




