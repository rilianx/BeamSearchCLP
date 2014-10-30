 #ifndef _BOX
 #define _BOX

#include "btBulletDynamicsCommon.h"
#include <BulletCollision/BroadphaseCollision/btDbvt.h>
#include <iostream>
#include <map>
#include <set>
#include <list>
#include "Policy.h"

using namespace std;
 namespace odp {
   #if defined(_MSC_VER) || defined(__BORLANDC__)
   typedef unsigned __int64 ulong64;
   typedef signed __int64 long64;
   #else
   typedef unsigned long long ulong64;
   typedef signed long long long64;
   #endif
   
   
   
  
  //Puede ser necesario optimizar metodo: void _removegap(btDbvtNode* node); (directamente desde nodo)
  class Box : public btDbvtAabbMm{
  
    public:
    //Nodo correspondiente
    btDbvtNode* node;
    long64 volume;

    
		Box() : btDbvtAabbMm(btDbvtAabbMm::FromMM (btVector3 (0,0,0), 
		btVector3 (0,0,0))), node(NULL){}
		
    Box(btDbvtNode* node) : btDbvtAabbMm(node->volume), node(node), volume(getW()*getL()*getH()) {}
	
    Box(btDbvtAabbMm volume, btDbvtNode* node) : btDbvtAabbMm(volume), node(node), volume(getW()*getL()*getH()) {}
    
	Box(long64 x1, long64 y1, long64 z1, long64 x2, long64 y2, long64 z2) : 
	    btDbvtAabbMm(btDbvtAabbMm::FromMM (btVector3 ((btScalar)x1,(btScalar)y1,(btScalar)z1), 
		btVector3 ((btScalar)x2,(btScalar)y2,(btScalar)z2))), volume(getW()*getL()*getH()), node(NULL){}
	
    bool operator>=(const Box &b) const {
      return (getW()>=b.getW() &&  getL()>=b.getL() && getH()>=b.getH());
    }

    bool operator<=(const Box &b) const {
      return (getW()<=b.getW() &&  getL()<=b.getL() && getH()<=b.getH());
    }

    bool includes(const Box &b) const {
      return (getX()<=b.getX() &&  getY()<=b.getY() && getZ()<=b.getZ() &&
       b.getXmax()<=getXmax() &&  b.getYmax()<=getYmax() && b.getZmax()<=getZmax() );
    }

 
    inline long64 getX() const{
      return (long64) (Mins().getX()+0.5);
    }

    inline long64 getY() const{
      return (long64) (Mins().getY()+0.5);
    }

    inline long64 getZ() const{
      return (long64) (Mins().getZ()+0.5);
    }

    inline long64 getXmax() const{
      return (long64) (Maxs().getX()+0.5);
    }

    inline long64 getYmax() const{
      return (long64) (Maxs().getY()+0.5);
    }

    inline long64 getZmax() const{
      return (long64) (Maxs().getZ()+0.5);
    }

    inline long64 getW() const{
      return (long64) (Lengths().getY()+0.5);
    }

    inline long64 getL() const{
      return (long64) (Lengths().getX()+0.5);
    }

    inline long64 getH() const{
      return (long64) (Lengths().getZ()+0.5);
    }

    inline btDbvtAabbMm get_inbox(){
        return btDbvtAabbMm::FromMM (Mins() + btVector3(0.5,0.5,0.5), Maxs() - btVector3(0.5,0.5,0.5));
    }
    
    inline btDbvtAabbMm get_outbox(){
        return btDbvtAabbMm::FromMM (Mins() - btVector3(0.5,0.5,0.5), Maxs() + btVector3(0.5,0.5,0.5));
    }
    
    inline long64 get_volume() const{
        return volume;
    }
    
    inline bool collides(const Box &b){
	  return !( getX()+getL()<=b.getX() || b.getX()+b.getL()<=getX() || getY()+getW()<=b.getY() || b.getY()+b.getW()<=getY() || getZ()+getH()<=b.getZ() || b.getZ()+b.getH()<=getZ() );
	}
   
    virtual void print() const;
    
  };
  

  
  //ordenar por anchor_dist (lexicografico), luego volumen
  class FSpace : public Box {
    public:
    

    list<long64> anchor_dist;
    bool anchor_corner[3]; 
	long64 manhattan_distance;
    
	FSpace(const FSpace& fsp, btDbvtNode* node) : Box(node), 
	anchor_dist(fsp.anchor_dist), manhattan_distance(fsp.manhattan_distance){
		anchor_corner[0]=fsp.anchor_corner[0]; 
		anchor_corner[1]=fsp.anchor_corner[1]; 
		anchor_corner[2]=fsp.anchor_corner[2];
	}
	
    FSpace(btDbvtAabbMm volume, btDbvtNode* node, long64 L, long64 W, long64 H, bool fsb=false) : Box(volume, node), manhattan_distance(0) {
        long64 dist;
        if(getX() <= L - (getX() + getL())){
		  dist=getX();
          anchor_corner[0]=0;			
        }else{
					dist= L - (getX() + getL());
          anchor_corner[0]=1;
        }
				anchor_dist.push_back(dist); 
				manhattan_distance+=dist;
       
       if(getY() <= W - (getY() + getW())){
				 dist=getY();
         anchor_corner[1]=0;
       }else{
				 dist=W - (getY() + getW());
         anchor_corner[1]=1;
       }
			anchor_dist.push_back(dist); 
			manhattan_distance+=dist;
       
       if(fsb || (getZ() <= H - (getZ() + getH()))){
			dist=getZ();
            anchor_corner[2]=0;
       }else{
			dist=H - (getZ() + getH());
            anchor_corner[2]=1;
       }
			anchor_dist.push_back(dist); 
			manhattan_distance+=dist;			 
			 

       anchor_dist.sort();
    }
    
    bool operator==(const FSpace& s)const{
       return(getX()==s.getX() && getY()==s.getY() && getZ()==s.getZ() && getW()==s.getW() && getL()==s.getL() && getH()==s.getH());
    }
    
    virtual void print() const{
      cout << "(FreeSpace)  Dim: " << getL() << " x " << getW() << " x " << getH() << "    ";
      cout << " Location: (" << getX() <<"," << getY() <<"," << getZ() <<")";
      list<long64> l=anchor_dist;
      list<long64> :: iterator it=l.begin();
      cout << " Anchor_dist:" << *it; it++;
      cout << "," << *it; it++;
      cout << "," << *it << "  volume:"<< get_volume() << endl;
    }
    
  };
  
  
 enum Orientation{LWH=0, LHW, WLH, WHL, HLW, HWL};

/* struct compareBoxTypes{
   bool operator() (BoxType* const& b1, BoxType* const& b2) const {
	   
	   
	if(b1->occupied_volume!=b2->occupied_volume)  return (b1->occupied_volume > b2->occupied_volume);
	// if(b1->volume!=b2->volume) return (b1->volume > b2->volume);
	 if(b1->l!=b2->l) return (b1->l>b2->l);
     if(b1->w!=b2->w) return (b1->w>b2->w);
     if(b1->h!=b2->h) return (b1->h>b2->h);
     if(less_than(b1->nb_boxes, b2->nb_boxes)) return true;
     else return false;
  }
 };
 */
 
 class BoxType{
 public:
	static int n;
    map<Orientation, bool> is_valid_orientation;
	int sort_dim[3];
	int id;
	int minl, minw, minh;
	
    BoxType(long64 l, long64 w, long64 h, int id, bool rot1=true, bool rot2=true, bool rot3=true);
   
    //get oriented dimensions of the box
    inline long64 getL(Orientation o){ 
      if(o==WLH || o==WHL)
	return w;
      else if(o==LWH || o==LHW)
	return l;      
      else if(o==HLW || o==HWL)
	return h;  
	  
	  return -1;      
    }
    
    inline long64 getW(Orientation o){ 
      if(o==LWH || o==HWL)
	return w;
      else if(o==WLH || o==HLW)
	return l;      
      else if(o==WHL || o==LHW)
	return h;    
	  	  return -1;     
    }
    
    inline long64 getH(Orientation o){ 
      if(o==HLW || o==LHW)
	return w;
      else if(o==HWL || o==WHL)
	return l;      
      else if(o==LWH || o==WLH)
	return h;   
	  	  return -1;            
    }
   
    bool operator<(const BoxType& bt) const{
    	return (id<bt.id);
    }
         
       long64 volume;

       void print();
        
       long64 w;
       long64 l;
       long64 h;   
   
 };
  

    void add(const Box& b1, const Box& b2, list<Box>& res);
	void subtract(const Box& b1, const Box& b2, list<Box>& res, long64 pack_l=-1, long64 pack_w=-1);
	long64 maxspace_vol(Box& b1, Box& b2);
    void filter_nonmaximalGaps(list<btDbvtAabbMm>& gaps);
	void filter_nonmaximalGaps(list<Box>& gaps);
    void print(const btDbvtAabbMm & b);
    
  
}

  
#endif
