#include "Block.h"

 using namespace odp;
 namespace odp{
 
 double Block::expon=2.0;
 long64 Block::calc=0;
 
 bool less_than(const map<BoxType*,int>& l, const map<BoxType*,int>& r)
 {
   // same types, proceed to compare maps here

   if(l.size() != r.size())
     return (l.size() < r.size());  // differing sizes, they are not the same

   map<BoxType*,int>::const_iterator i, j;
   for(i = l.begin(), j = r.begin(); i != l.end(); ++i, ++j)
   { 
     if(i->first != j->first)
       return (*(i->first) < *(j->first));
	 if(i->second!=j->second)
	   return (i->second < j->second);
   }

   return false; //they are equal
 }



//Crea una lista de a lo mas 3 bloques juntando dos bloques (sin rotacion)
list<Block*> Block::createBlocks(Block& b1, Block& b2, long64 maxL, long64 maxW, 
long64 maxH, double min_fr, bool fsb){
	
	list<Block*> blocks;
	
    long64 volume=b1.volume+b2.volume;
	long64 occupied_volume=b1.occupied_volume+b2.occupied_volume;
	long64 vol=0;
			
	map<BoxType*,int>* nb_boxes=NULL;
	
		
	map<BoxType*,int>::const_iterator it_nb;	
	
  for(int i=0; i<3; i++){

	long64 ll= max(b1.l,b2.l); 
	long64 ww= max(b1.w,b2.w); 
	long64 hh= max(b1.h,b2.h);
    
    long64 pack_area_l;
    long64 pack_area_w;
    
	switch(i){
		case 0: 
		  if(fsb && (b1.h!=b2.h || b1.pack_area_l!=b1.l || b2.pack_area_l!=b2.l)) continue;
		  pack_area_l=b1.pack_area_l+b2.pack_area_l;
		  pack_area_w=min(b1.pack_area_w,b2.pack_area_w);
		  ll=(b1.l+b2.l); break;
		case 1: 
		  if(fsb && (b1.h!=b2.h || b1.pack_area_w!=b1.w || b2.pack_area_w!=b2.w)) continue;
		  pack_area_l=min(b1.pack_area_l,b2.pack_area_l);
		  pack_area_w=b1.pack_area_w+b2.pack_area_w;
		  ww=(b1.w+b2.w); break;
		case 2: 
		  if(fsb && ((b1.pack_area_w<b2.w || b1.pack_area_l<b2.l) && 
		  (b2.pack_area_w<b1.w || b2.pack_area_l<b1.l))) continue;
		  pack_area_l=min(b1.pack_area_l,b2.pack_area_l);
		  pack_area_w=min(b1.pack_area_w,b2.pack_area_w);		  
		  hh=(b1.h+b2.h);
	}
 
	if(ll>maxL || ww>maxW || hh>maxH) continue;
	
	long64 vol= ll*ww*hh;
	if((double) occupied_volume / (double) vol >= min_fr ){
			
		if(!nb_boxes){
			nb_boxes=new map<BoxType*,int>();
            for(it_nb = b1.nb_boxes.begin(); it_nb!=b1.nb_boxes.end(); it_nb++)
            (*nb_boxes)[(*it_nb).first] += (*it_nb).second;
			
            for(it_nb = b2.nb_boxes.begin(); it_nb!=b2.nb_boxes.end(); it_nb++)
            (*nb_boxes)[(*it_nb).first] += (*it_nb).second;		
	    }
			
		Block* new_block=new Block(ll, ww, hh, b1.occupied_volume + b2.occupied_volume, *nb_boxes);
		new_block->pack_area_w=pack_area_w;
		new_block->pack_area_l=pack_area_l;
		
		//se ponen uno sobre otro
		if(fsb && i==2){
			if(b1.pack_area_w >= b2.w && b1.pack_area_l >= b2.l){
			  new_block->block1=&b1;  //bloque inferior
			  new_block->block2=&b2;  //bloque superior
		    }else{
			  new_block->block1=&b2;  //bloque inferior
			  new_block->block2=&b1;  //bloque superior
		    }
		}else{
			  new_block->block1=&b1; 
			  new_block->block2=&b2; 			
		}


	    //~ //ubicacion relativa del segundo bloque
	    switch(i){
			case 0: new_block->x2=new_block->block1->l; break;
			case 1: new_block->y2=new_block->block1->w; break;
			case 2: new_block->z2=new_block->block1->h;
	    }

		blocks.push_back(new_block);						
        
	}
 }
 if(nb_boxes) delete nb_boxes;
 return blocks;
	 
}




void get_feasible_blocks(list<Block*>& blocks, map<BoxType*,int>& nb_left_boxes, list<Block*>& feasible_blocks){
   list<Block*>:: iterator it=blocks.begin();
   for(;it!=blocks.end();it++)
      if((*it)->feasible(nb_left_boxes)) feasible_blocks.push_back(*it);

}


    
}
