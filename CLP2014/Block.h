 #ifndef _BLOCK
 #define _BLOCK

#include "btBulletDynamicsCommon.h"
#include <BulletCollision/BroadphaseCollision/btDbvt.h>
#include <iostream>
#include <map>
#include <set>
#include <list>
#include "Policy.h"
#include "Box.h"
#include <fstream>
#include <vector>


using namespace std;
 namespace odp {



	 bool less_than(const map<BoxType*,int>& l, const map<BoxType*,int>& r);

	 //Importante: no admiten rotación
  class Block{
    public:
    
    int id;
    long64 w,l,h;
    long64 occupied_volume;
    long64 volume;
    long64 sqr_volumes;
	int n_boxes;
	
    //tipo de cajas, número de cajas
    map<BoxType*,int> nb_boxes;
    long64 fitness; 
    long64 ranking;   

	 
	//para impresion de resultados...
	Block* block1; //ubicación relativa al block contenedor: (0,0,0)
	Block* block2; long64 x2,y2,z2;
	 
	//usados en variante con restriccion fsb (full supported blocks)
	long64 pack_area_w;
	long64 pack_area_l;
	 
    static double expon;
   
    
   
    //only for set search
    Block(FSpace& fspace) : w(fspace.getW()), l(fspace.getL()), h(fspace.getH()), volume(fspace.volume)
    , block1(NULL), block2(NULL), x2(0),y2(0),z2(0), pack_area_w(0), pack_area_l(0) {
    	
    }
    
    static long64 calc;
    
    void print_loc(long64 x0, long64 y0,long64 z0,ofstream& myfile, vector<Box*>& boxes){
		if(n_boxes > 1){
			if(block1 ==NULL) cout << "error" << endl;
			block1->print_loc(x0, y0, z0,myfile, boxes);
			block2->print_loc(x0+x2, y0+y2, z0+z2,myfile,boxes);
			
			//~ cout << "block " << id << ": (" << x0 << "," << y0 << "," << z0 << "); " << 
			 //~ ": (" << x0+l << "," << y0+w << "," << z0+h << "); " << endl;
		}else{
			myfile << nb_boxes.begin()->first->id; 
			myfile << "	(" << x0 << "," << y0 << "," << z0 << ")	" << 
			 "(" << x0+l << "," << y0+w << "," << z0+h << ")" << endl;
			boxes.push_back(new Box(x0, y0, z0, x0+l, y0+w, z0+h));
			 calc+=l*w*h;
		}
		
		
	}
	
    //Simple Block (nW x nL x nH)
    Block(BoxType* bt, Orientation o=LWH, int nL=1, int nW=1, int nH=1) :  l(nL*bt->getL(o)), w(nW*bt->getW(o)), 
    h(nH*bt->getH(o)), volume(l*w*h),block1(NULL), block2(NULL), x2(0),y2(0),z2(0), 
    pack_area_w(nW*w), pack_area_l(nL*l){
         nb_boxes[bt]=nL*nW*nH;
		 n_boxes=nb_boxes[bt];
         occupied_volume=l*w*h; //caso particular 100% ocupado (simple block)
		 sqr_volumes=(long64) (n_boxes*pow((double) bt->volume,expon));
		 
    }
    
	Block(int l, int w, int h, long64 occupied_volume, map<BoxType*,int>& nb_boxes) : l(l), w(w), h(h), 
	occupied_volume(occupied_volume), nb_boxes(nb_boxes), volume(l*w*h),   
	block1(NULL),  block2(NULL), x2(0),y2(0),z2(0), pack_area_w(0), pack_area_l(0){

	   n_boxes=0; sqr_volumes=0;
       map<BoxType*,int>::iterator it_nb;
       for(it_nb = nb_boxes.begin(); it_nb!=nb_boxes.end(); it_nb++){
       	  n_boxes+=(*it_nb).second;
		  sqr_volumes+=(long64) ((*it_nb).second*pow((*it_nb).first->volume,expon));
	   }
       
   }
	
   bool operator<=(const Box &b) const {
       return (w<=b.getW() &&  l<=b.getL() && h<=b.getH());
   }
	
   //Crea una lista de a lo mas 3 bloques juntando dos bloques (sin rotacion)
   //fsb=full supported blocks
   static list<Block*> createBlocks(Block& b1, Block& b2, long64 maxL, long64 maxW, long64 maxH, double min_fr=0.98, 
   bool fsb=false);
    
    //Block(const Block &b) : l(b.l), w(b.w), h(b.h), occupied_volume(b.occupied_volume), volume(b.volume),
    //nb_boxes(b.nb_boxes), fitness(b.fitness), id(b.id) {}
    
    
    Box get_location(const FSpace &free_space){ 
       long64 minx=free_space.getX(), miny=free_space.getY(), minz=free_space.getZ();
       
       /**** coloca el bloque en anchor_corner *****/
       if(free_space.anchor_corner[0]==1) minx=free_space.getX()+free_space.getL() - l;
       if(free_space.anchor_corner[1]==1) miny=free_space.getY()+free_space.getW() - w;
       if(free_space.anchor_corner[2]==1) minz=free_space.getZ()+free_space.getH() - h;
       
       btDbvtAabbMm vol= btDbvtAabbMm::FromMM (btVector3((btScalar) minx,(btScalar) miny,(btScalar) minz),
          btVector3 ((btScalar) (minx+l),(btScalar) (miny+w),(btScalar) (minz+h)));
       return Box(vol, NULL);
    }
    
    //retorna si se puede crear el bloque con las cajas que quedan
    bool feasible(map<BoxType*,int>& nb_left_boxes){
       map<BoxType*,int>::iterator it_nb;
       for(it_nb = nb_boxes.begin(); it_nb!=nb_boxes.end(); it_nb++)
         if(nb_left_boxes[(*it_nb).first] < (*it_nb).second) return false;
       
       return true;
    }
    
    void print() const{
      cout << l << " " << w << " " << h << "  volume:" << occupied_volume 
		  << "("<< ((double) occupied_volume/(double)volume)<<")" <<  endl;
       map<BoxType*,int>::const_iterator it_nb;
       for(it_nb = nb_boxes.begin(); it_nb!=nb_boxes.end(); it_nb++){
       cout << "   * " << it_nb->first->id << " " <<  it_nb->first->l << " " <<  it_nb->first->w << " "<<
        it_nb->first->h << ":" << it_nb->second  << "  volume:" << (it_nb->first->volume*it_nb->second) << endl;
       }
    }  

    void eval_volume_fitness(){
        fitness = volume;
    }
    

  };
  
  struct compareBlocks{
    bool operator() (Block* const& b1, Block* const& b2) const {
		if(b1->occupied_volume!=b2->occupied_volume)  return (b1->occupied_volume > b2->occupied_volume);
		// if(b1->volume!=b2->volume) return (b1->volume > b2->volume);
		 if(b1->l!=b2->l) return (b1->l>b2->l);
	     if(b1->w!=b2->w) return (b1->w>b2->w);
	     if(b1->h!=b2->h) return (b1->h>b2->h);
	     if(less_than(b1->nb_boxes, b2->nb_boxes)) return true;
	     else return false;
   }
  };
  

  
  void get_feasible_blocks(list<Block*>& blocks, map<BoxType*,int>& nb_left_boxes, list<Block*>& feasible_blocks);

  struct by_fitness{

    //funcion usada para ordenar bloques
    bool operator() (Block* const& b1, Block* const& b2) const{
      if(b1->fitness != b2->fitness) return (b1->fitness > b2->fitness );
      else return b1->id<b2->id;
    }
  };
  
  
  struct by_id{

    //funcion usada para ordenar bloques
    bool operator() (Block* const& b1, Block* const& b2) const{
      return b1->id<b2->id;
    }
  };
	
  struct by_ranking{

    //funcion usada para ordenar bloques
    bool operator() (Block* const& b1, Block* const& b2) const{
      if(b1->ranking != b2->ranking) return (b1->ranking > b2->ranking );
      else return b1<b2;
    }
  };
	

  //funcion usada para ordenar los tipos de cajas
  inline bool by_volume(Block* & b1, Block* & b2){
	  return (b1->occupied_volume>b2->occupied_volume);
  }
  
  inline bool by_volume_div_nb_box(Block* & b1, Block* & b2){
	  int vdn1=int((double)b1->occupied_volume/(double)b1->n_boxes+0.99999);
	  int vdn2=int((double)b2->occupied_volume/(double)b2->n_boxes+0.99999);
	  if(vdn1!=vdn2) return (vdn1>vdn2);
	  	  
	  return (b1->occupied_volume>b2->occupied_volume);
  } 
  
  inline bool by_nb_box(Block* & b1, Block* & b2){	 
	  if(b1->n_boxes != b2->n_boxes) 	return (b1->n_boxes<b2->n_boxes);  
	  return (b1->occupied_volume>b2->occupied_volume);
  }
  
  
  
  
  //funcion usada para ordenar los tipos de cajas
  inline bool by_sqr_volumes(Block* & b1, Block* & b2){
	  return (b1->sqr_volumes>b2->sqr_volumes);
  }


  
  //objeto bloque // no tipo
  class Block_obj : public Box {
	public:
	
	
	Block* block; //tipo de bloque
 	//cantidad de bloques que soporta
	//bloques que no soportan otros bloques pueden ser descartados
	int n_supported_boxes;
	//bloques de soporte, estos bloques no pueden ser descartados mientras el este colocado el bloque actual
	Block_obj* support_block[3]; 
	bool anchor_corner[3]; 
	
	Block_obj() {};
	
	Block_obj(Block_obj* b) : Box(b->node), block(b->block), n_supported_boxes(0) {
		support_block[0]=NULL;
		support_block[1]=NULL;
		support_block[2]=NULL;	
		anchor_corner[0]=b->anchor_corner[0];
		anchor_corner[1]=b->anchor_corner[1];
		anchor_corner[2]=b->anchor_corner[2];
		node=NULL;
	}
	
    Block_obj(Block* block, Block_obj* block2) : 
    Box(
        (block2->anchor_corner[0]==1)? block2->getX()+block2->getL() - block->l:block2->getX(),
        (block2->anchor_corner[1]==1)? block2->getY()+block2->getW() - block->w:block2->getY(),
        (block2->anchor_corner[2]==1)? block2->getZ()+block2->getH() - block->h:block2->getZ(),
        (block2->anchor_corner[0]==1)? block2->getX()+block2->getL():block2->getX() + block->l,
        (block2->anchor_corner[1]==1)? block2->getY()+block2->getW():block2->getY() + block->w,
        (block2->anchor_corner[2]==1)? block2->getZ()+block2->getH():block2->getZ() + block->h), block(block), n_supported_boxes(0){ 
		anchor_corner[0]=block2->anchor_corner[0];
		anchor_corner[1]=block2->anchor_corner[1];
		anchor_corner[2]=block2->anchor_corner[2];
			
 		support_block[0]=NULL;
		support_block[1]=NULL;
		support_block[2]=NULL;      
        node=NULL;
    }

	  
    Block_obj(Block* block, const FSpace &free_space) : 
    Box(
        (free_space.anchor_corner[0]==1)? free_space.getX()+free_space.getL() - block->l:free_space.getX(),
        (free_space.anchor_corner[1]==1)? free_space.getY()+free_space.getW() - block->w:free_space.getY(),
        (free_space.anchor_corner[2]==1)? free_space.getZ()+free_space.getH() - block->h:free_space.getZ(),
        (free_space.anchor_corner[0]==1)? free_space.getX()+free_space.getL():free_space.getX() + block->l,
        (free_space.anchor_corner[1]==1)? free_space.getY()+free_space.getW():free_space.getY() + block->w,
        (free_space.anchor_corner[2]==1)? free_space.getZ()+free_space.getH():free_space.getZ() + block->h), block(block), n_supported_boxes(0){ 
		anchor_corner[0]=free_space.anchor_corner[0];
		anchor_corner[1]=free_space.anchor_corner[1];
		anchor_corner[2]=free_space.anchor_corner[2];
			
 		support_block[0]=NULL;
		support_block[1]=NULL;
		support_block[2]=NULL;     
		node=NULL; 
       
    }
    
  };


}

#endif
