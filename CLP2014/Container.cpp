
#include "Container.h"
#include "Box.h"
#include "btBulletDynamicsCommon.h"
#include <BulletCollision/BroadphaseCollision/btDbvt.h>
#include <iostream>
#include <map>
#include <set>

namespace odp{

    double Container::alpha=1.0;
    bool Container::fsb=false;
    
    Container::Container(int L, int W, int H) : W(W), H(H), L(L), _occupied_volume0(0) { 
      btDbvtAabbMm first_gap=btDbvtAabbMm::FromMM (btVector3 (0.0,0.0,0.0), btVector3 ((btScalar) L,(btScalar) W, (btScalar) H));
      btDbvtNode * node=_tree_fspaces.insert(first_gap, NULL);
      node2space[node]=free_spaces.insert(FSpace(node->volume,node,L,W,H)).first;
      mL= new int[L+1];
      mH= new int[H+1];
      mW= new int[W+1];
    }
    

    Container::Container(Container& cont): L(cont.L), W(cont.W), H(cont.H), mL(cont.mL), mW(cont.mW), mH(cont.mH),
     _occupied_volume0(cont.occupied_volume()){
	    
	    //se copian los free spaces
        set<FSpace, by_manhattan_distance>::iterator it=cont.free_spaces.begin();
	    while(it!=cont.free_spaces.end()){
            btDbvtNode * node=_tree_fspaces.insert(*it, NULL);
            node2space[node]=free_spaces.insert(FSpace(*it,node)).first;
	        it++;
	    }
	
        //se copian los bloques (locations)
	    list<Block_obj*>::iterator it2=cont.locations.begin();
	    map<Block_obj*,Block_obj*> created;

	    while(it2!=cont.locations.end()){
			if(!created[*it2]){
			
  			   Block_obj* blox =  create_bloxR(*it2,created);
               if(blox->n_supported_boxes==0)
                 removable_blocks.insert(blox);
			}
            it2++;
	    }
        _occupied_volume0=0;
    }

    Block_obj* Container::create_bloxR(Block_obj* b, map<Block_obj*,Block_obj*> &created){
		Block_obj* blox =  new Block_obj(*b);
		created[b]=blox;
		
		btDbvtNode * node=_tree_blocks.insert(*blox, NULL);
		node2block[node]=blox;
		locations.push_back(blox);
		//_occupied_volume0-=blox->block->occupied_volume;
		
		blox->node=node;
	    for(int i=0;i<3;i++)
	       if(b->support_block[i]){
			  if(created[b->support_block[i]]) blox->support_block[i]=created[b->support_block[i]];
	          else blox->support_block[i]=create_bloxR(b->support_block[i],created);
		   }
	    return blox;
    }

    void Container::print_fspaces(){
	set<FSpace, by_manhattan_distance> :: iterator it = free_spaces.begin();
	
        for(; it!=free_spaces.end();it++){
           (*it).print();
	}
    }
    
    void Container::print_blocks(){
      list<Block_obj*> :: iterator it=locations.begin();
//       cout << gaps.size() << endl;
      for(;it!=locations.end();it++)
	(*it)->block->print();
      
    }
    
    bool by_volume2(Box & b1, Box & b2) {
     return (b1.volume>b2.volume) ;/* || (v1==v2 && compare()(Aabb_to_point(b1),Aabb_to_point(b2)));*/
    }
	
//     struct volume_sort{
     bool volume_sort(const btDbvtAabbMm& b1, const btDbvtAabbMm& b2) {
	btScalar v1=b1.Lengths().getX()*b1.Lengths().getY()*b1.Lengths().getZ();
	btScalar v2=b2.Lengths().getX()*b2.Lengths().getY()*b2.Lengths().getZ();
	
      return v1<v2 ;/* || (v1==v2 && compare()(Aabb_to_point(b1),Aabb_to_point(b2)));*/
    }
//     };
         

    void Container::_remove_fspace(btDbvtNode* node){ 
	    free_spaces.erase(FSpace(node->volume,node,L,W,H));
	   _tree_fspaces.remove(node);
    }

    bool Container::is_feasible(const Box& fsp, map<BoxType*,int>& nb_left_boxes){
        map<BoxType*,int>::iterator bl=nb_left_boxes.begin();
        for(;bl!=nb_left_boxes.end();bl++){
           if(bl->second==0) continue;
           if(bl->first->volume > fsp.volume) continue;
             
           for(int i=0;i<6;i++){
             if(bl->first->is_valid_orientation[Orientation(i)])
               if(bl->first->getL(Orientation(i)) <= fsp.getL() && bl->first->getW(Orientation(i)) <= fsp.getW() 
               && bl->first->getH(Orientation(i)) <= fsp.getH() ){
                 return true;
               }
           }
        }		
		return false;
    }

    long64 Container::remove_unfeasible_spaces(map<BoxType*,int>& nb_left_boxes){
	   long64 total_volume=0;
       set<FSpace, by_manhattan_distance>::iterator fs=free_spaces.begin();
       while(free_spaces.end()!=fs){
          if(!is_feasible(*fs,nb_left_boxes)){
              //set<FSpace, by_manhattan_distance>::iterator fs2=fs; fs2++;
    		   node2space.erase(fs->node);
			  _tree_fspaces.remove(fs->node);
			  free_spaces.erase(fs++);
              //fs=fs2;
          }else{
			 total_volume+=fs->volume;
			 fs++;
		  }
      }
	  return total_volume;
    }   

  long64 Container::compute_VLoss(Box& box, Box& free_space){
      int tmp1, tmp2;
    
      int resFront = mL[((tmp1=int(box.getX()) - int(free_space.getX())) >0)? tmp1:0];
      int resBack  = mL[((tmp2=int(free_space.getXmax()) - int(box.getXmax())) >0)? tmp2:0];
      if(tmp1<0) tmp1=0; if(tmp2<0) tmp2=0;
      int bl = (free_space.getL() - tmp1 - tmp2);
      
      int resLeft = mW[((tmp1=int(box.getY()) - int(free_space.getY())) >0)? tmp1:0];
      int resRight= mW[((tmp2=int(free_space.getYmax()) - int(box.getYmax())) >0)? tmp2:0];
      if(tmp1<0) tmp1=0; if(tmp2<0) tmp2=0;
      int bw = (free_space.getW() - tmp1 - tmp2);

      int resTop    = mH[((tmp1=int(box.getZ()) - int(free_space.getZ())) >0)? tmp1:0];
      int resBottom = mH[((tmp2=int(free_space.getZmax()) - int(box.getZmax())) >0)? tmp2:0];
      if(tmp1<0) tmp1=0; if(tmp2<0) tmp2=0;   
      int bh = (free_space.getH() - tmp1 - tmp2);
 
      int realL=min(mL[resFront] + mL[resBack] + bl,mL[free_space.getL()]);
      int realW=min(mW[resLeft] + mW[resRight] + bw,mW[free_space.getW()]);
      int realH=min(mH[resTop] + mH[resBottom] + bh,mH[free_space.getH()]);
      
      return (mL[free_space.getL()]*mW[free_space.getW()]*mH[free_space.getH()] -  long64(realL)*long64(realW)*long64(realH));
   }

   long64 Container::compute_VLossCorrected(Block* block, const FSpace& free_space){
      int resL=mL[int(free_space.getL())] - int(block->l);
      int resW=mW[int(free_space.getW())] - int(block->w);
      int resH=mH[int(free_space.getH())] - int(block->h);
      
      if(resL<0 || resW<0 || resH<0) return 0;
      
      int xl = resL - mL[resL];
      int xw = resW - mW[resW];
      int xh = resH - mH[resH];

      return block->occupied_volume - (block->l*xh*xw + block->w*xl*xh + block->h*xw*xl + 
        block->l*block->w*xh + block->l*xw*block->h + xl*block->w*block->h + xl*xw*xh);
   }

   long64 Container::compute_VLossFitness(Block* block, const FSpace& free_space){
      long64 resL=free_space.getL() - block->l;
      long64 resW=free_space.getW() - block->w;
      long64 resH=free_space.getH() - block->h;
      
      if(resL<0 || resW<0 || resH<0) return -1e10;
      
      int realL=block->l + mL[resL];
      int realW=block->w + mW[resW];
      int realH=block->h + mH[resH];

      return block->occupied_volume - 
      alpha*(mL[free_space.getL()]*mW[free_space.getW()]*mH[free_space.getH()] - long64(realL)*long64(realW)*long64(realH));
      
      //return 1e14 + block->occupied_volume - alpha*(free_space.volume - long64(realL)*long64(realW)*long64(realH));
   }
   
   void Container::solveKnapsack(map<BoxType*,int>& nb_left_boxes){
		

    map<BoxType*,int>::iterator it_nb;
    int maxDim=max(L,max(W,H));
		
    bool flag[maxDim+1]; flag[0]=true;

	
    for(int i=1; i<=L; i++) flag[i]=false;	
    set<int> ws;
    for(it_nb = nb_left_boxes.begin(); it_nb!=nb_left_boxes.end(); it_nb++){
	int n=(*it_nb).second;
			
        if(n > 0){
                
             for(int ii=0; ii<3; ii++){
                 int w=L+1;
                 switch(ii){
                       case 0: w=int((*it_nb).first->l); break;
                       case 1: 
                         if((*it_nb).first->is_valid_orientation[WLH] || (*it_nb).first->is_valid_orientation[WHL]) 
                           w=int((*it_nb).first->w);
                         break;
                       case 2:
                         if((*it_nb).first->is_valid_orientation[HLW] || (*it_nb).first->is_valid_orientation[HWL]) 
                           w=int((*it_nb).first->h);
                         break;        
                       default:  break;               
                 }
                   
                 if(w>L || ws.find(w)!=ws.end()) continue;                  
                 ws.insert(w);
			
		 for(int i=0; i<=L-w; i++){
                     if(!flag[i]) continue;
                   
                     for(int j=min(n,(int(L)-i)/w); j>=1; j--)
                        flag[i+ w*j]=true;
                 }

            }
			  
       }
   }
   
   for(int i=0; i<=L; i++){
      if(flag[i]) mL[i]=i;
      else mL[i]=mL[i-1];
      //~ cout << "mL[" << i << "]:" << mL[i] << endl;
   }


    for(int i=1; i<=W; i++) flag[i]=false;	
    ws.clear();
    for(it_nb = nb_left_boxes.begin(); it_nb!=nb_left_boxes.end(); it_nb++){
	int n=(*it_nb).second;
			
        if(n > 0){
                
             for(int ii=0; ii<3; ii++){
                 int w=W+1;
                 switch(ii){
                       case 0: w=int((*it_nb).first->w); break;
                       case 1: 
                         if((*it_nb).first->is_valid_orientation[HLW] || (*it_nb).first->is_valid_orientation[WLH]) 
                           w=int((*it_nb).first->l);
                         break;
                       case 2:
                         if((*it_nb).first->is_valid_orientation[LHW] || (*it_nb).first->is_valid_orientation[WHL]) 
                           w=int((*it_nb).first->h);
                         break;        
                       default:  break;               
                 }
                   
                 if(w>W || ws.find(w)!=ws.end()) continue;                  
                 ws.insert(w);
			
		 for(int i=0; i<=W-w; i++){
                     if(!flag[i]) continue;
 
                                  
                     for(int j=min(n,(int(W)-i)/w); j>=1; j--)
                        flag[i+ w*j]=true;
                 }

            }
			  
       }
   }
   
   for(int i=0; i<=W; i++){
      if(flag[i]) mW[i]=i;
      else mW[i]=mW[i-1];
			//~ cout << "mW[" << i << "]:" << mW[i] << endl;
      
   }
   
    for(int i=1; i<=H; i++) flag[i]=false;	
    ws.clear();
    for(it_nb = nb_left_boxes.begin(); it_nb!=nb_left_boxes.end(); it_nb++){
	int n=(*it_nb).second;
			
        if(n > 0){
                
             for(int ii=0; ii<3; ii++){
                 int w=H+1;
                 switch(ii){
                       case 0: w=int((*it_nb).first->h); break;
                       case 1: 
                         if((*it_nb).first->is_valid_orientation[WHL] || (*it_nb).first->is_valid_orientation[HWL]) 
                           w=int((*it_nb).first->l);
                         break;
                       case 2:
                         if((*it_nb).first->is_valid_orientation[HLW] || (*it_nb).first->is_valid_orientation[LHW]) 
                           w=int((*it_nb).first->w);
                         break;        
                       default:  break;               
                 }

                 if(w>H || ws.find(w)!=ws.end()) continue;                  
                 ws.insert(w);
			
		 for(int i=0; i<=H-w; i++){
                     if(!flag[i]) continue;
                   
                     for(int j=min(n,(int(H)-i)/w); j>=1; j--)
                        flag[i+ w*j]=true;
                 }

            }
			  
       }
   }
   

   for(int i=0; i<=H; i++){
      if(flag[i]) mH[i]=i;
      else mH[i]=mH[i-1];
      //~ cout << "mH[" << i << "]:" << mH[i] << endl;
   }
	 //~ exit(0);
   
}

	

	void Container::block1(btDbvtAabbMm& _outbox, btDbvtAabbMm& _inbox, Box& blox, 
	list<Box>& new_gaps, Block* block){
        //Get gaps colliding with the box
        Policy p(Policy::COLLISION_SET);
        _tree_fspaces.collideTV ( _tree_fspaces.m_root, _outbox, p);
		
	  
        while(!p.collided_nodes.empty()){
  	       btDbvtNode* node=p.collided_nodes.front();
  	       p.collided_nodes.pop_front();
  	       Box b(node->volume,NULL);	   
  	       if(Intersect(node->volume,_inbox)){
			   if(block) subtract(b, blox, new_gaps, block->pack_area_l, block->pack_area_w);
			   else subtract(b, blox, new_gaps);
  	      }else
  	         new_gaps.push_back(b);
		   
 
   		   free_spaces.erase(node2space[node]);
   		   node2space.erase(node);
   		   _tree_fspaces.remove(node);
        }			
	}
	
	//elimina un bloque
    void Container::remove(Block_obj* blox, map<BoxType*,int>& nb_left_boxes){	
	
	  Block* block=blox->block;
	  btDbvtAabbMm _outbox= blox->get_outbox(); 	

	  locations.remove(blox);

      //Se actualiza la cantidad de cajas restantes
      map<BoxType*,int>::iterator it_nb;

      for(it_nb = block->nb_boxes.begin(); it_nb!=block->nb_boxes.end(); it_nb++)
        nb_left_boxes[(*it_nb).first]+=(*it_nb).second;

      
      //ACTUALIZACION de FREE SPACES
      //se recuperan los agujeros colindantes con el bloque a eliminar
      //se generan nuevos agujeros maximales, estos son almacenados en outputs
      //se eliminan los agujeros colindantes
      //se agregan los agujeros generados
      list<Box> outputs;
      outputs.push_back(*blox);

      //Get gaps colliding with the box (inputs)
      Policy p(Policy::COLLISION_SET);
      _tree_fspaces.collideTV ( _tree_fspaces.m_root, _outbox, p);
	          
      while(!p.collided_nodes.empty()){
           
	       btDbvtNode* s2=p.collided_nodes.front();

           list<Box>::iterator tmp=outputs.begin();
 	       insert_spaces(outputs,*node2space[s2],tmp); 

  	       for(list<Box>::iterator s1=outputs.begin();s1!=outputs.end();){

			   list<Box> nb;
			   add(*s1,*node2space[s2],nb);
			   
               list<Box>::iterator ref=s1;
			   bool removed_ref=insert_spaces(outputs,nb,ref); 
			   
			   if(!removed_ref) s1++;
			   else s1=ref;
		   }

  	       p.collided_nodes.pop_front();
  	       //se elimina el espacio
 		   free_spaces.erase(node2space[s2]);
   		   node2space.erase(s2);
   		   _tree_fspaces.remove(s2);	

	    }

    //cout << "outputs:" <<  outputs.size() << endl;
  	for(list<Box>::iterator s1=outputs.begin();s1!=outputs.end();s1++){
	 	btDbvtNode * node = _tree_fspaces.insert(*s1, NULL);
	    node2space[node]=free_spaces.insert(FSpace(node->volume,node,L,W,H,fsb)).first;		
    }
		

  	//se actualizan soportes 	 
      for (int i=0;i<3; i++){
		  Block_obj* b=blox->support_block[i];
		  blox->support_block[i]=NULL;
		  if(b){
		    b->n_supported_boxes--;
		    if(b->n_supported_boxes==0) removable_blocks.insert(b);
		  }
	  }

	  removable_blocks.erase(blox);
   	  node2block.erase(blox->node);
   	  _tree_blocks.remove(blox->node);        

     
	  
   
	}


// si el nuevo espacio (new_box) está contenido en algún espacio de initial, no hace nada
// si el nuevo espacio no está contenido en ningún espacio de initial, se agrega a initial
//    todos los espacios contenidos por new_box son eliminados de initial
    bool Container::insert_spaces(list<Box>& initial, const Box& new_box,list<Box>::iterator &ref){
			bool flag = false;
			bool removed_ref=false;
		    for(list<Box>::iterator s=initial.begin();s!=initial.end();){
				if(s->includes(new_box)) {flag=true;  break;}
				if(new_box.includes(*s)) {
				   if(s==ref){
				      s=initial.erase(s);	
				      ref=s; 
				      removed_ref=true;
				   }else
				      s=initial.erase(s);	
				   	
	               //break;
				}else
				  s++;
			}
			
			if (!flag)
	               initial.push_back(new_box);				
	    return removed_ref;
	}

    void Container::insert_spaces(list<Box>& initial, list<Box>& news, list<Box>::iterator &ref){
		bool removed_ref=false;
		for(list<Box>::iterator s2=news.begin();s2!=news.end();s2++)
			removed_ref |= insert_spaces(initial, *s2, ref);
		return removed_ref
	}
    



    //Coloca el bloque en el espacio libre indicado
    int Container::insert(Block_obj* blox, map<BoxType*,int>& nb_left_boxes, bool remspaces){
		
      //Block_obj* blox=new Block_obj(blo);
      Block* block=blox->block;
      btDbvtAabbMm _inbox= blox->get_inbox();
      btDbvtAabbMm _outbox= blox->get_outbox();    
      
      locations.push_front(blox);
      //blocks.push_front(block);
      
      //Se actualiza la cantidad de cajas restantes
      map<BoxType*,int>::iterator it_nb;
      for(it_nb = block->nb_boxes.begin(); it_nb!=block->nb_boxes.end(); it_nb++)
        nb_left_boxes[(*it_nb).first]-=(*it_nb).second;

      list<Box> new_gaps;

      //se elimina el bloque (blox) de la estructura de free spaces
      block1(_outbox, _inbox, *blox, new_gaps, (fsb)? block:NULL);

 
      new_gaps.sort(by_volume2);
      filter_nonmaximalGaps(new_gaps);

      //if(save) cout << new_gaps.size() << endl;
      //add the resultant gaps
      while(!new_gaps.empty()){	
		if(!remspaces || is_feasible(new_gaps.front(),nb_left_boxes)){
	       btDbvtNode * node = _tree_fspaces.insert(new_gaps.front(), NULL);
	       node2space[node]=free_spaces.insert(FSpace(node->volume,node,L,W,H,fsb)).first;
	   }
	   new_gaps.pop_front();
      }
      

      //Se buscan los bloques de soporte del nuevo bloque y se actualizan
      Policy p(Policy::COLLISION_SET);
      _tree_blocks.collideTV (_tree_blocks.m_root, _outbox, p);     

      while(!p.collided_nodes.empty()){
  	       btDbvtNode* nn=p.collided_nodes.front(); //se pueden ordenar por tamanno?
  	       Block_obj* b=node2block[nn];
  	       int i=-1;
  	       if(!(blox->getX()==0 || blox->getXmax()==L) && !blox->support_block[0] && (b->getX() == blox->getXmax() || b->getXmax() == blox->getX()) )
               i=0;
		   else if(!(blox->getY()==0 || blox->getYmax()==W) && !blox->support_block[1] && (b->getY() == blox->getYmax() || b->getYmax() == blox->getY()) )
		       i=1;
		   else if(!(blox->getZ()==0 || blox->getZmax()==H) && !blox->support_block[2] && (b->getZ() == blox->getZmax() || b->getZmax() == blox->getZ()) )
		       i=2;
           
           if(i!=-1){
		       blox->support_block[i]=b;
			   if(b->n_supported_boxes==0) removable_blocks.erase(b);
			   b->n_supported_boxes++;
		   }
  	       
  	       p.collided_nodes.pop_front();
	  }
    
      btDbvtNode * node = _tree_blocks.insert(*blox, NULL);
      node2block[node]=blox;
      blox->node=node;
      removable_blocks.insert(blox);
      
      return OK;
      
    }
    
    void initialize_rankings(list<Block*>& blocks, map<BoxType*,int>& nb_left_boxes, set<FSpace, by_manhattan_distance>& fspaces){
 	   list<Block*>::iterator bl=blocks.begin();
 	   for(;blocks.end()!=bl; bl++){
 	      (*bl)->ranking=0;
 	      bool ok=false;

 	      if(!(*bl)->feasible(nb_left_boxes)){
 		    	(*bl)->ranking=-2e10;
 		       continue;
 		    }

 	      set<FSpace, by_manhattan_distance>::iterator fs=fspaces.begin();
	      
 	      for(;fs!=fspaces.end() && !ok;fs++){
 	         if((*bl)->volume > fs->volume) continue;
               
 		     if(**bl<=*fs) {ok=true; break;}
 	      }
            
 	      if(!ok){
 		    (*bl)->ranking=-2e10;
 	      }
 	  }	    
    }
    
    void remove_unfeasible_blocks_rank(list<Block*>& blocks){
	    list<Block*>::iterator bl=blocks.begin();
	    for(;blocks.end()!=bl;){
		    if((*bl)->ranking==-2e10){
		    	bl=blocks.erase(bl);
			continue;
		    }
		    bl++;
	    }
    	
    }
 
 	void Container::remove_unfeasible_blocks_conserv(list<Block*>& blocks, map<BoxType*,int>& nb_left_boxes, set<FSpace, by_manhattan_distance>& fspaces){
	   list<Block*>::iterator bl=blocks.begin();
	   while(blocks.end()!=bl){
	      (*bl)->ranking=0;

	      if(!(*bl)->feasible(nb_left_boxes)){
			(*bl)->ranking=-2e10;
	            bl=blocks.erase(bl);
		      continue;
		  }
		
         bl++;

	  }

    }
    
	void Container::remove_unfeasible_blocks(list<Block*>& blocks, map<BoxType*,int>& nb_left_boxes, set<FSpace, by_manhattan_distance>& fspaces){
	   list<Block*>::iterator bl=blocks.begin();
	   while(blocks.end()!=bl){
	      (*bl)->ranking=0;
	      bool ok=false;

	      if(!(*bl)->feasible(nb_left_boxes)){
			(*bl)->ranking=-2e10;
	            bl=blocks.erase(bl);
		      continue;
		}
		
          
	      set<FSpace, by_manhattan_distance>::iterator fs=fspaces.begin();
	  
	      for(;fs!=fspaces.end() && !ok;fs++){
	         if((*bl)->volume > fs->volume) continue;
             
		     if(**bl<=*fs) {ok=true; break;}
	          
	      }

	      if(!ok){
		    (*bl)->ranking=-2e10;
	          bl=blocks.erase(bl);
	      }else
			 bl++;

	  }

      }
      
    bool Container::fit_spaces(Block_obj& b1, Block_obj& b2){
		
	     
	     Block_obj b12(b1.block,&b2); //b1 in b2 space
	     Block_obj b21(b2.block,&b1); //b2 in b1 space
	     
	     btDbvtAabbMm _inbox1= b12.get_inbox();  
	     btDbvtAabbMm _inbox2= b21.get_inbox(); 
	      
	     Policy p(Policy::COLLISION_SET);
         _tree_blocks.collideTV (_tree_blocks.m_root, _inbox1, p);
         
         if(p.collided_nodes.size()==1){
			 p.collided_nodes.clear();
	         _tree_blocks.collideTV (_tree_blocks.m_root, _inbox2, p);
	         if(p.collided_nodes.size()==1) return true;	 
		 }
		 return false;
         
	}

	
	void remove_unfeasible_blocks(list<Block*>& blocks, map<BoxType*,int>& nb_left_boxes){
	   list<Block*>:: iterator it=blocks.begin();
	   while(it!=blocks.end()){
	      if(!(*it)->feasible(nb_left_boxes)){
	         it=blocks.erase(it);
	      }else it++;
	   }
	}
    
}
  
