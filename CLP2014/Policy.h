#ifndef _POLICY
#define _POLICY

#include "btBulletDynamicsCommon.h"
#include <BulletCollision/BroadphaseCollision/btDbvt.h>
// #include "Box.h"
#include <list>

using namespace std;
namespace odp {

  class is_included: public exception
{
  virtual const char* what() const throw()
  {
    return "My exception happened";
  }
};
  
      bool Contain(btDbvtAabbMm& b, btDbvtAabbMm& a);
class Policy : public btDbvt::ICollide {

  
    
public:
  enum type{COUNT_ONLY, COLLISION_SET, INCLUDED_SET, ALL};
  int nb_collisions;
  type t; 
  
  list<btDbvtNode*> collided_nodes;
  list<btDbvtNode*> is_included_in_nodes;
  list<btDbvtNode*> include_nodes;
  
  btDbvtAabbMm* referential_box;

  
  Policy(type t) : t(t), nb_collisions(0), referential_box(NULL) {}
  virtual void 	Process (const btDbvtNode *, const btDbvtNode *){}
  virtual void 	Process (const btDbvtNode *n){
    nb_collisions++;
    if(t==COLLISION_SET)
      collided_nodes.push_back(const_cast<btDbvtNode*>(n));
    else if(t==INCLUDED_SET){
      if(referential_box->Contain(n->volume)) collided_nodes.push_back(const_cast<btDbvtNode*>(n));
    }else if(t==ALL){
      if(Contain(const_cast<btDbvtNode*>(n)->volume,*referential_box)) {
	is_included_in_nodes.push_back(const_cast<btDbvtNode*>(n));
	throw is_included();
      }
      else if(Contain(*referential_box,const_cast<btDbvtNode*>(n)->volume)) include_nodes.push_back(const_cast<btDbvtNode*>(n));   
      else collided_nodes.push_back(const_cast<btDbvtNode*>(n));
    }
     
  }
  virtual void 	Process (const btDbvtNode *n, btScalar s){}
  virtual bool 	Descent (const btDbvtNode *){return false;}
  virtual bool 	AllLeaves (const btDbvtNode *){return false;}
  void clear(){
    nb_collisions=0;
    collided_nodes.clear();
    is_included_in_nodes.clear();
    include_nodes.clear();
//     if(referential_box) delete referential_box;
//     referential_box=NULL;
  }
  
  ~Policy(){
//       if(referential_box) delete referential_box;
  }
  
};

}

#endif