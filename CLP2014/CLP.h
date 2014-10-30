 #ifndef _CLP
 #define _CLP

#include "btBulletDynamicsCommon.h"
#include <BulletCollision/BroadphaseCollision/btDbvt.h>
#include <iostream>
#include <map>
#include <set>
#include <list>
#include "Policy.h"
#include "Block.h"

using namespace std;
 namespace odp {
   #if defined(_MSC_VER) || defined(__BORLANDC__)
   typedef unsigned __int64 ulong64;
   typedef signed __int64 long64;
   #else
   typedef unsigned long long ulong64;
   typedef signed long long long64;
   #endif
   
   //container loading problem
   class CLP {
      public:
     
      ~CLP(){
         nb_boxes.clear();

         while(box_types.size()!=0){

            delete box_types.front();
            box_types.pop_front();
         }
      }
	  
      CLP(ifstream& in, double total_time=100);
     
      void add_box_type(BoxType* btype, int nb){
         box_types.push_back(btype);
         nb_boxes[btype]=nb;
      }
     
	  set<Block*, compareBlocks> single_box_block_generator();
	  set<Block*, compareBlocks> general_block_generator(double min_fr=1.0, int max_bl=10000, bool fsb=false);
	 
	  inline double get_time(){
		  return (double(clock()-begin_time)/double(CLOCKS_PER_SEC));
	  }
	  
	  inline bool timeout(){
	  	return (get_time()>=total_time);
	  }
	 
      long64 W;
      long64 L;
      long64 H;
	  double total_time;
	  clock_t begin_time;
     
      list<BoxType*> box_types;
      map<BoxType*,int> nb_boxes;
   };
   
   
   
   
}

#endif
