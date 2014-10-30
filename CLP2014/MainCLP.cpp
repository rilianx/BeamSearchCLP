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
		
		Node(list<Action>& path, long64 est_vol, long64 vol, Container* cont, map<BoxType*,int>* nb_left_boxes, Node* father) : 
		path(path), est_vol(est_vol), 
		vol(vol), cont(cont), nb_left_boxes(nb_left_boxes), father(father){
			
		}
		
		//root node
		Node(Container* cont, map<BoxType*,int>* nb_left_boxes) : cont(cont), 
		nb_left_boxes(nb_left_boxes), est_vol(0), vol(0), depth(0), father(NULL){
			
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
    long64 est_vol; //volume estimado (ub) 
    long64 vol; //volumen ocupado
		int depth;
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
	
	 vertex& operator*( const double& d){
 		vertex* ret=new vertex(size);
             for( int i = 0; i < size; i++ )
                 ret->x[i] = d * x[i];
 		return *ret;	 	
	 }
	
	vertex& operator+( const vertex& vv  ){
		vertex* ret=new vertex(vv.size);
            for( int i = 0; i < size; i++ )
                ret->x[i] = x[i] + vv.x[i];
		return *ret;		
	}

	vertex& operator-( const vertex& vv  ){
		vertex* ret=new vertex(vv.size);
            for( int i = 0; i < size; i++ )
                ret->x[i] = x[i] - vv[i];
		return *ret;		
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
	
	void print() {
		cout << "(";
		for(int i=0; i<size; i++)
		   cout << x[i] << ((i+1==size)? "":",");
		cout <<")" << endl;
		
	}
	
};

void _2LAw(CLP& clp, set<Block*, compareBlocks>& set_blocks, vector<SortedList*>& sorted_lists, 
    vector<vertex>& weights, int n, long64 &best_volume, double grasp, list<Node>& best_path);
void update_path_GRASP(CLP& clp, list<Block*>& blocks, vector<SortedList*>& sorted_lists, vector<vertex>& weights,
Container& cont, map<BoxType*,int>& nb_left_boxes, int w, long64& best_volume, list<Node>& best_path);

long64 update_path_kLA(CLP& clp, list<Block*>& blocks, vector<SortedList*>& sorted_lists, Node& node, vector<vertex>& weights, int w, list<Node*>* childs, list<Action>& path, int depth, int k);

long64 greedy(list<Block*>& blocks, Container& cont, map<BoxType*,int>& nb_left_boxes, 
	double grasp_perc=1.0, list<Action>* path=NULL);
Block* select_block_2LA(CLP& clp, list<Block*>& blocks, Container& cont, map<BoxType*,int>& nb_left_boxes, int w, long64& max_vol);
long64 greedy(
  list<Block*>& blocks, vector<SortedList*>& sorted_lists, double *w, Container& cont, map<BoxType*,int>& nb_left_boxes, list<Action>* path, double grasp=0.0);
vertex nelder_mead(CLP& clp, set<Block*, compareBlocks> set_blocks, vector<SortedList*>& sorted_lists, int s, double grasp, list<Node>& best_path,
long64& best);
void random_grasp(CLP& clp, vector<SortedList*>& sorted_lists, double grasp, list<Action>& best_path);
void beam_search(CLP& clp, set<Block*, compareBlocks>& set_blocks, vector<SortedList*>& sorted_lists, 
    vector<vertex>& weights, int w, long64 &best_volume, double grasp, list<Action>& best_path, bool parallel_search=false);


int greedy_calls=0;
double GR_PERC=0.9;
int gruns=0;



int main(int argc, char** argv)
{


	ifstream in(argv[1]);
	int _inst0=atoi(argv[2]); //instancia de partida (e.g. 25)
	int _inst1=atoi(argv[3]); //instancia de partida (e.g. 25)
	long64 w0=atoi(argv[4]);
  srand(atoi(argv[5]));
	bool parallel=(atoi(argv[6])==1);


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
	   CLP clp(in,20); //se inicializa el problema, 600 segundos de tiempo
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
		     weights, w, best_volume, 0, best_path, parallel);
			//	cout << w << ":"<< double(best_volume) / double(clp.W*clp.L*clp.H) << " " << (clp.get_time()-time0) << endl;
				w=ceil(double(w)*sqrt(2));	

	   }
		    	      
     //cout << double(best_volume) << endl;
		 cout << w <<":" << double(best_volume) / double(clp.W*clp.L*clp.H) << " " << (clp.get_time()-time0) << endl;
         while(set_blocks.size()>0){
   		     delete *set_blocks.begin();
   		     set_blocks.erase(set_blocks.begin());
         }
	   
      }




}


inline bool by_eval(vertex & v1, vertex & v2){
  return (v1.eval>v2.eval);
}


double run_greedy(CLP& clp, set<Block*, compareBlocks>& set_blocks, vector<SortedList*>& sorted_lists, double* weights, int n, long64 &best, 
double grasp, list<Action>& best_path){
	//gruns+=n;
	long64 total_volume=0;
	long64 best_vol=0;
	for(int i=0; i<n; i++){
	   Container cont(clp.L,clp.W,clp.H);
	   map<BoxType*,int> nb_left_boxes=clp.nb_boxes;
         list<Block*> blocks(set_blocks.begin(), set_blocks.end());
         initialize_rankings(blocks, nb_left_boxes, cont.free_spaces);	
	
         list<Block*>::iterator bl=blocks.begin();
         for(int i=0;bl!=blocks.end();i++, bl++) (*bl)->id=i; 	
	   for(int i=0; i<sorted_lists.size(); i++)  sorted_lists[i]->init(blocks);
	   cont.solveKnapsack(nb_left_boxes);

         list<Action> new_path;	
	   long64 volume=greedy(blocks, sorted_lists, weights, cont, nb_left_boxes, &new_path, grasp);

	   if(volume>best){
		   best_path=new_path;
		   best=volume;
	   } 
	   
	   if(volume>best_vol)
		   best_vol=volume;
	 //  cout << double(volume) / double(clp.W*clp.L*clp.H) << endl;
	   total_volume+=volume;
     }
     return (double(total_volume) / (double) n) / double(clp.W*clp.L*clp.H);
}

void random_grasp(CLP& clp, vector<SortedList*>& sorted_lists, double grasp, list<Action>& best_path){
	int n=gruns;
	 set<Block*, compareBlocks> set_blocks = clp.general_block_generator();	
	vertex x0(sorted_lists.size());
	long64 best=0;
	
	for(int i=0; i<n; i++){
	    for(int j=0; j<sorted_lists.size(); j++)
      	        x0[j]=double(rand())/double(RAND_MAX);
	    run_greedy(clp, set_blocks, sorted_lists, x0.x, 1, best, grasp, best_path);
      }
	cout << "grasp:" << double(best) / double(clp.W*clp.L*clp.H) << endl;
}

//Algoritmo nelder_mead se usa para encontrar estrategia relativamente buena
vertex nelder_mead(CLP& clp, set<Block*, compareBlocks> set_blocks, vector<SortedList*>& sorted_lists, int s, double grasp, list<Action>& best_path,
long64& best){
      int n=sorted_lists.size();

	double dn=double(n);
	double a=1.0;
	
	vertex x0(n);
	for(int j=0; j<n; j++)
      	x0[j]=double(rand())/double(RAND_MAX); //A number between 0 and 1

	x0.eval=run_greedy(clp, set_blocks, sorted_lists, x0.x, s, best, grasp, best_path);

	list<vertex> vertices;
	
	vertices.push_back(x0);
	
	
	double p = a/(dn*sqrt(2))*(sqrt(dn+1)+dn-1);
	double q = a/(dn*sqrt(2))*(sqrt(dn+1)-1);
	cout << "p:" << p << endl;
	cout << "q:" << q << endl;
	
	vertex xi(n);
	for(int i=1; i<n+1; i++){
		for(int j=0; j<n; j++)
			xi[j]=x0[j]+((i-1==j)? p:0) + ((i-1 != j)? q:0);
		
		
		xi.eval=run_greedy(clp, set_blocks, sorted_lists, xi.x, s, best, grasp, best_path);
		vertices.push_back(xi);
	}

      //NELDER-MEAD
      double alpha=1, beta=0.5, gamma=2.0, delta=0.5;
      for(int iter=0; iter<10; iter++){	
	  vertices.sort(by_eval);
	  list<vertex>:: iterator it=vertices.begin();
	//    for(;it!=vertices.end(); it++)
	  //	  cout << it->eval << ",";
	   // cout << endl;
	  //Calculo del centroide
	  vertex xc(n);
	  it=vertices.begin();
	  for(int i=0; i<vertices.size()-1; it++, i++)	  xc = xc + (*it);
	  for(int j=0; j<n; j++) xc[j]/=double(n);	

	  vertex& xl=*vertices.begin();
	  it=vertices.end();
	  it--; vertex& xh=*it;
	  it--; vertex& xs=*it;

        if(xl.eval - xh.eval <0.001) break;


       // xl.print();

        vertex xr=xc+(xc-xh)*alpha;
	  xr.eval=run_greedy(clp, set_blocks, sorted_lists, xr.x, s, best, grasp, best_path);
	  //cout << xr.eval << endl;
	  //Expansion
        if(xr.eval > xl.eval){ //xr es el mejor
		 // cout << "expansion" << endl;
      	vertex xe= xc+(xr-xc)*gamma;
		xe.eval=run_greedy(clp, set_blocks, sorted_lists, xe.x, s, best, grasp, best_path);
		if(xe.eval>xr.eval) 
			xh=xe;
		else
			xh=xr;	
	  }else if(xr.eval > xs.eval){ //mejor que el segundo peor
		 // cout << "simply accepted" << endl;
		xh=xr;
	  }else{//Contraction
		  //cout << "contraction" << endl;
	    if(xr.eval > xh.eval) xh=xr;
	    vertex xcc=xc+(xh-xc)*beta;
	    xcc.eval=run_greedy(clp, set_blocks, sorted_lists, xcc.x, s, best, grasp, best_path);
	   // cout << xcc.eval << endl;
	    if(xcc.eval<xh.eval){ //Shrink
		    //cout << "shrink" << endl;
		    
		it=vertices.begin();
	    	for(it++; it!=vertices.end(); it++){
			vertex &v=*it;
			v = (v + xl)*0.5;
	  	      v.eval=run_greedy(clp, set_blocks, sorted_lists, v.x, s, best, grasp, best_path);
	    	}	    
	    }else{
		    xh=xcc;
	    }
		
        }
      }
	
	cout << "Nelder_mead:" << double(best) / double(clp.W*clp.L*clp.H) << endl;
      return *vertices.begin();
	
	
	//cout << "best:" << double(best) / double(clp.W*clp.L*clp.H) << endl;
}

bool by_est_vol(Node*  n1, Node*  n2){
		if(n1->est_vol != n2->est_vol) return (n1->est_vol > n2->est_vol );
		return (n1->vol < n2->vol);
    
  }


void beam_search(CLP& clp, set<Block*, compareBlocks>& set_blocks, vector<SortedList*>& sorted_lists, 
    vector<vertex>& weights, int w, long64 &best_volume, double grasp, list<Action>& best_path, bool parallel_search){
			int node_limit=w;
			
	    Container* cont= new Container(clp.L,clp.W,clp.H);
	    map<BoxType*,int>* nb_left_boxes= new map<BoxType*,int>(clp.nb_boxes);
	    list<Block*> blocks(set_blocks.begin(), set_blocks.end());
			for(int i=0; i<sorted_lists.size(); i++)  sorted_lists[i]->init(blocks);    
			
			list<Node*> L; list<Node*> newL;
			list<Node*> allNodes; //para segunda etapa->mejora
			L.push_back(new Node(cont, nb_left_boxes));
			
			bool first=true;
			int depth=1;
			while(L.size()!=0){
				
				list<Node*>::iterator itL=L.begin();
				
				for(;itL!=L.end();itL++){
				  
			   	Node* node=*itL; 
					list<Action> best_path;
						long64 volume;
					if(!parallel_search){
						volume=update_path_kLA(clp, blocks, sorted_lists, *node, weights, w, &newL, best_path, 0, 0);
					}else{
						list<Node*> childL;
						volume=update_path_kLA(clp, blocks, sorted_lists, *node, weights, w, &childL, best_path, 0, 0);
						childL.sort(by_est_vol); 
						 if(childL.size()>0)
							     newL.push_back(childL.front());
						
						list<Node*>::iterator it=childL.begin();		
            if(childL.size()>0) it++;
						 for(;it!=childL.end();it++)
							 	 	delete *it;		
					}
					
						 
					if(volume>best_volume)
						best_volume=volume;
			  }
				
			  first=false;
				 if(!parallel_search){//eliminar identicos y los peores

				//cout << newL.size() << endl;

					newL.sort(by_est_vol); 
					list<Node*>::iterator it=newL.begin();		
					long64 	old_est_vol, old_vol;	
					int i=0;  
					for(;it!=newL.end();it++, i++){
						if(i>=node_limit) {delete *it; continue;}				
						if(i>0 && (*it)->est_vol==old_est_vol && (*it)->vol==old_vol){ delete *it; it=newL.erase(it);	it--; i--; continue;}
						old_est_vol=(*it)->est_vol;
						old_vol=(*it)->vol;
						 
					}					 
					if(newL.size()>node_limit) newL.resize(node_limit);
				}
				
				
				
				//cout << "children:" << endl;
				 list<Node*>::iterator it=newL.begin();	
				for(;it!=newL.end();it++){
					//cout << double((*it)->est_vol) / double(clp.W*clp.L*clp.H)<< endl;
					  //se copian cont y nb_left_boxes del nodo padre y se modifican
	 				  Container* cont_tmp= new Container(*((*it)->cont));
	 				  map<BoxType*,int>* boxes_tmp=new map<BoxType*,int>(*((*it)->nb_left_boxes));
							
					    cont_tmp->insert((*it)->action.block, (*it)->action.x, (*it)->action.y, 
					 	   (*it)->action.z, *boxes_tmp, false);	
							
						//	cont_tmp->remove_unfeasible_spaces(*boxes_tmp);
				
											
							
						
					  (*it)->instantiate(cont_tmp, boxes_tmp);
						(*it)->depth=depth;
				}
									
				for(it=L.begin();it!=L.end();it++){
					if((*it)->cont) delete (*it)->cont;
					if((*it)->nb_left_boxes) delete (*it)->nb_left_boxes;
					allNodes.push_back(*it);			
					//delete *it;
				}
		      
				
		
				L=newL;
				newL.clear();
				depth++;
			}
			
			list<Node*>::iterator it;
			for(it=allNodes.begin();it!=allNodes.end();it++){
		//		cout << (*it)->depth << ":" << (*it)->vol << " " << (*it)->est_vol << endl;
			   delete *it;
		  }
			
		}






   //retorna numero de bloques que calzan en el espacio correspondiente
   int actualize_rankings_fspace(list<Block*>& blocks, const FSpace& fspace, map<BoxType*,int>& nb_left_boxes){
   	int nb_blocks=0;
      list<Block*>::iterator it=blocks.begin();
   	for(;it!=blocks.end();it++){
   	      if((*it)->ranking==-2.0) continue;	
            if(!(*it)->feasible(nb_left_boxes)){
                   (*it)->ranking=-2.0;
      		 continue;		    
   		}else if(!(**it<=fspace))
   		    (*it)->ranking=-1.0;
   		else{
   	          (*it)->ranking=0.0;
   		    nb_blocks++;
   	      }		
      }
      return nb_blocks; 	
   }

   void selectBlocks_rank(list<Block*>& blocks, vector<SortedList*>& sorted_lists, double *weights, Container& cont, 
   const FSpace& fspace, int n, map<BoxType*,int>& nb_left_boxes, set<Block*, by_ranking>& blocks_tmp){
	      
	blocks_tmp.clear();   
	
   	list<Block*>:: iterator it=blocks.begin();
   	for(;it!=blocks.end();it++){
	   if((*it)->ranking < 0) continue;

   	   if(blocks_tmp.size()<n || (*it)->ranking  < (*blocks_tmp.rbegin())->ranking){

     	      if(blocks_tmp.size()==n){
       	     set<Block*, by_ranking>::iterator ii=blocks_tmp.end(); --ii;
     	           blocks_tmp.erase(ii);
     	      }
   	      blocks_tmp.insert(*it);

   	   }
      }

     	      
   }

   Block* selectBlock_greedyrank(list<Block*>& blocks, int nb_blocks, vector<SortedList*>& sorted_lists, 
   double *w, Container& cont, FSpace& fspace,  map<BoxType*,int>& nb_left_boxes, double perc){
	   
	  
	  int size=max(int(perc*double(nb_blocks)+0.9999),1);
	  int r=rand()%size;

        set<Block*, by_ranking> blocks_tmp;
	  selectBlocks_rank(blocks, sorted_lists, w, cont, 
	      fspace, r+1, nb_left_boxes, blocks_tmp);

	  set<Block*, by_ranking>::iterator it=blocks_tmp.end();
	  it--;
	  return *(it);
	  	  
   }
   
   Block* selectBlock_rank(list<Block*>& blocks, int nb_blocks, vector<SortedList*>& sorted_lists, 
   double *w, Container& cont, FSpace& fspace,  map<BoxType*,int>& nb_left_boxes){

	   return selectBlock_greedyrank(blocks, nb_blocks, sorted_lists, 
	      w, cont, fspace,  nb_left_boxes, 0);
   }
   

   
  
  long64 greedy(
	  list<Block*>& blocks, vector<SortedList*>& sorted_lists, double *w, Container& cont, 
        map<BoxType*,int>& nb_left_boxes, list<Action>* path, double grasp){
	      
         while(!cont.free_spaces.empty()){
	
       	   FSpace fspace = *(cont.free_spaces.begin()); //K3: min anchor corn		   
		  
     	        //se coloca en el fspace la primera caja que quepa   
					 int nb_blocks=actualize_rankings_fspace(blocks, fspace, nb_left_boxes);
					 if(nb_blocks!=0){
						 actualize_rank(nb_blocks, blocks, sorted_lists, w, nb_left_boxes, cont, &fspace);
						 Block* block=selectBlock_greedyrank(blocks, nb_blocks, sorted_lists, w, cont, fspace, nb_left_boxes, grasp);

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


     //retorna mejor ruta (eventualmente podria retornar lista 
     //de mejores rutas) y mejor volumen encontrado
  long64 update_path_kLA(CLP& clp, list<Block*>& blocks, vector<SortedList*>& sorted_lists, Node& node, 
	vector<vertex>& weights, int w, list<Node*>* childs, list<Action>& path, int depth, int k){
		//cout << depth << endl;
		//se selecciona el primer espacio "factible" según anchor_distance
		
		if(depth==0) node.cont->solveKnapsack(*node.nb_left_boxes);
		
 	  set<FSpace, by_anchor_distance>::iterator itsp=node.cont->free_spaces.begin();	
 	  int nb_blocks=0;     
    while(nb_blocks==0 && itsp!=node.cont->free_spaces.end()){
			
 		     initialize_rankings(blocks, *node.nb_left_boxes, node.cont->free_spaces);
 		     nb_blocks=actualize_rankings_fspace(blocks, *itsp, *node.nb_left_boxes);
		     
 		     if(nb_blocks==0) {
			     
					 if(depth==0){
							node.cont->remove_fspace(*itsp);
							itsp=node.cont->free_spaces.begin();	
			 		 }else
					    itsp++;	 
		     }		           
 	  }
		if(nb_blocks==0) return 0;

		FSpace fspace = (*itsp);
		
		set<Block*, by_id> blocks_tmp; //bloques en el nivel depth del árbol
		
    //se seleccionan hasta w bloques con cada heurística de orden
		for(int i=0; i<weights.size(); i++){
			set<Block*, by_ranking> blocks_tmpp;
     	          //calcula ranking de las cajas restantes
			actualize_rank(nb_blocks, blocks, sorted_lists, weights[i].x,       
				*node.nb_left_boxes, *node.cont, &fspace);
			selectBlocks_rank(blocks, sorted_lists, weights[i].x, *node.cont, fspace, w, *node.nb_left_boxes, blocks_tmpp);
			blocks_tmp.insert(blocks_tmpp.begin(), blocks_tmpp.end());
							
		}
	/*	
		if(depth==0 && node.path.size()>0){
			Block* b=node.path.front().block;
			if(!(*b <= fspace)){
				cout << " error!"  << endl;	
			 }
	    else blocks_tmp.insert(b);
		}
		*/
						

		long64 best_volume=node.cont->occupied_volume();
		//long64 best_volume=node.est_volume;

		set<Block*, by_id>:: iterator it=blocks_tmp.begin();

		for(int i=-1; it!=blocks_tmp.end() && !clp.timeout(); i++){
	    
			if(i==-1 && (depth!=0 || node.path.size()==0)) i=0;
			
      Block* b;
			
			if(i==-1){
				//follow the best_path
				b=node.path.front().block;
		    node.cont->insert(b, node.path.front().x, node.path.front().y, 
		 	   node.path.front().z, *node.nb_left_boxes);		
				
			}else if(i==-2){
  		   node.cont->free_spaces.erase(node.cont->node2space[fspace.node]);
  		   node.cont->node2space.erase(fspace.node);
  		   node.cont->_tree_fspaces.remove(fspace.node);
				//node.cont->remove_fspace(fspae);
				
			}else{
				b=*it;
			  node.cont->insert(b, fspace, *node.nb_left_boxes);
			}
	    //for(int j=0; j<weights.size(); j++){
	     //for(l=0..tries)
			list<Action> new_path;		    
			long64 est_volume=0;
			if(depth<k){
				est_volume=update_path_kLA(clp, blocks, sorted_lists, node, weights, 
					w, NULL, new_path, depth+1, k);
			}else if(k>=0){
				Container cont_tmp(*node.cont);
				map<BoxType*,int> boxes_tmp=*node.nb_left_boxes;
				initialize_rankings(blocks, boxes_tmp, cont_tmp.free_spaces);
				est_volume=greedy(blocks, sorted_lists, weights[0].x, cont_tmp, boxes_tmp, &new_path, 0);
		   		    	 
			}	//else k==-1 --> sólo generación de nodos




			Action action;
			if(i==-2 && est_volume>node.cont->occupied_volume()){ 
					//otherwise the node is final
				  action=new_path.front();
				  new_path.pop_front();
			}else if(i!=-2){
				action=Action(b, node.cont->locations.front().getX(),
						node.cont->locations.front().getY(), node.cont->locations.front().getZ());
			
				if(depth==0){
						
					if(i==-1 && node.est_vol > est_volume){ 
					   //se recupera mejor ruta a través del nodo
						new_path=node.path;
						new_path.pop_front();
						est_volume=node.est_vol;
					}

					Node* n= new Node(new_path, est_volume, node.cont->occupied_volume(), node.cont, node.nb_left_boxes, &node);
					n->action=action;
		
					childs->push_back(n);
				
				}

				if(est_volume > best_volume){
					best_volume=est_volume;
					new_path.push_front(action);
					path=new_path;
				}
			}

			  if(i!=-2)
			    node.cont->remove_last_block(*node.nb_left_boxes);
				else {
					node.cont->insert_fspace(fspace);
				}


				if(i>=0) it++;
		}


		return best_volume;
						
	}
  


	
  




