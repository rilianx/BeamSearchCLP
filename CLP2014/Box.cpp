#include "Box.h"
 
 using namespace odp;
 namespace odp{
   	 
    
	 
    BoxType::BoxType(long64 l, long64 w, long64 h, int id, bool rot1, bool rot2, bool rot3) : 
    l(l), w(w), h(h), volume(w*l*h), id(id){

//       cout << w << l << h << endl;
      if(!rot3) {
	    cout << "rot3 is false!" << endl;
	    exit(0);
      }
    
      is_valid_orientation[WLH]=rot3;
      is_valid_orientation[LWH]=rot3;

	  is_valid_orientation[WHL]=rot1;
	  is_valid_orientation[HWL]=rot1;

	  is_valid_orientation[HLW]=rot2;
	  is_valid_orientation[LHW]=rot2;
	  
	  
	  minl=l;
	  if((rot3||rot1) && w<minl) minl=w;
	  if((rot1||rot2) && h<minl) minl=h;
	  
	  minw=w;
	  if((rot3||rot2) && l<minw) minw=l;
	  if((rot1||rot2) && h<minw) minw=h;
	  
	  minh=h;
	  if((rot1) && l<minh) minh=l;
	  if((rot2) && w<minh) minh=h;
	  
    }
    
    
   

    void Box::print() const{
      cout << "  Dim: " << getL() << " x " << getW() << " x " << getH() << "    ";
      if(node)
       cout << "Location: (" << getX() <<"," << getY() <<"," << getZ() <<")";
      
      cout << endl;
    }
    


    long64 maxspace_vol(Box& b1, Box& b2){
        long64 b1_xmin= b1.getX();
        long64 b1_ymin= b1.getY();
        long64 b1_zmin= b1.getZ();
        long64 b2_xmin= b2.getX();
        long64 b2_ymin= b2.getY();
        long64 b2_zmin= b2.getZ();	
        long64 b1_xmax= b1.getX() + b1.getL();
        long64 b1_ymax= b1.getY() + b1.getW();
        long64 b1_zmax= b1.getZ() + b1.getH();
        long64 b2_xmax= b2.getX() + b2.getL();
        long64 b2_ymax= b2.getY() + b2.getW();
        long64 b2_zmax= b2.getZ() + b2.getH();
   
        long64 maxvol=0.0;
		
        if(b1_xmax-b2_xmax > b2_xmin-b1_xmin){
		   if(b1_xmax-b2_xmax > 0){
        	   long64 vol=(b1_xmax-b2_xmax)*b1.getW()*b1.getH();
			   maxvol=vol;
		   }
		}else{
		   if(b2_xmin-b1_xmin>0){
        	   long64 vol=(b2_xmin-b1_xmin)*b1.getW()*b1.getH();
			   maxvol=vol;
		   }
        }

        if(b1_ymax-b2_ymax>b2_ymin-b1_ymin){
			if(b1_ymax-b2_ymax>0){
        	    long64 vol=b1.getL()*(b1_ymax-b2_ymax)*b1.getH();
			    if(maxvol<vol) maxvol=vol;
			}
		}else{
			if(b2_ymin-b1_ymin>0){
        	    long64 vol=b1.getL()*(b2_ymin-b1_ymin)*b1.getH();
			    if(maxvol<vol) maxvol=vol;
		    }
        }
		
        if(b1_zmax-b2_zmax>b2_zmin-b1_zmin){
			if(b1_zmax-b2_zmax>0){
        	     long64 vol=b1.getL()*b1.getW()*(b1_zmax-b2_zmax);
			     if(maxvol<vol) maxvol=vol;
			}
		}else{
			if(b2_zmin-b1_zmin>0){
        	    long64 vol=b1.getL()*b1.getW()*(b2_zmin-b1_zmin);
			    if(maxvol<vol) maxvol=vol;
			}
        }
		
		return maxvol;

    }


   void add(const Box& b1, const Box& b2, list<Box>& res){
      

      long64 b1_xmin= b1.getX();
      long64 b1_ymin= b1.getY();
      long64 b1_zmin= b1.getZ();
      long64 b2_xmin= b2.getX();
      long64 b2_ymin= b2.getY();
      long64 b2_zmin= b2.getZ();	
      long64 b1_xmax= b1.getX() + b1.getL();
      long64 b1_ymax= b1.getY() + b1.getW();
      long64 b1_zmax= b1.getZ() + b1.getH();
      long64 b2_xmax= b2.getX() + b2.getL();
      long64 b2_ymax= b2.getY() + b2.getW();
      long64 b2_zmax= b2.getZ() + b2.getH();

      
      //interbox:
      long64 g_xmin = max(b1_xmin, b2_xmin);
      long64 g_xmax = min(b1_xmax, b2_xmax);     
      long64 g_ymin = max(b1_ymin, b2_ymin);
      long64 g_ymax = min(b1_ymax, b2_ymax);  
      long64 g_zmin = max(b1_zmin, b2_zmin);
      long64 g_zmax = min(b1_zmax, b2_zmax);
     
      if(g_xmax -g_xmin < 0 || g_ymax -g_ymin < 0 || g_zmax -g_zmin < 0 ) return; //no intersection
      
//       //the sum is equivalent to one of the boxes:
//       if(g_xmin==b1_xmin && g_xmax==b1_xmax && g_ymin==b1_ymin && g_ymax == b1_ymax && g_zmin== b1_zmin && g_zmax==b1_zmax) return;
//       if(g_xmin==b2_xmin && g_xmax==b2_xmax && g_ymin==b2_ymin && g_ymax == b2_ymax && g_zmin== b2_zmin && g_zmax==b2_zmax) return;
      
      long64 gexp_xmin= min(b1_xmin, b2_xmin);
      long64 gexp_xmax= max(b1_xmax, b2_xmax);      
      long64 gexp_ymin= min(b1_ymin, b2_ymin);
      long64 gexp_ymax= max(b1_ymax, b2_ymax);    
      long64 gexp_zmin= min(b1_zmin, b2_zmin);
      long64 gexp_zmax= max(b1_zmax, b2_zmax);   

      if(((b1_xmin<b2_xmin && b1_xmax < b2_xmax) || (b1_xmin>b2_xmin && b1_xmax > b2_xmax)) &&
	gexp_xmax-gexp_xmin>0 &&  g_ymax-g_ymin >0 &&  g_zmax-g_zmin >0 ){
	
	  res.push_back(Box(gexp_xmin,g_ymin,g_zmin,gexp_xmax,g_ymax,g_zmax));
      }
      

  
      if(((b1_ymin<b2_ymin && b1_ymax < b2_ymax) || (b1_ymin>b2_ymin && b1_ymax > b2_ymax)) &&
	g_xmax-g_xmin>0 &&  gexp_ymax-gexp_ymin >0 &&  g_zmax-g_zmin >0 ){

	  res.push_back(Box(g_xmin,gexp_ymin,g_zmin,g_xmax,gexp_ymax,g_zmax));
      }

       
      if(((b1_zmin<b2_zmin && b1_zmax < b2_zmax) || (b1_zmin>b2_zmin && b1_zmax > b2_zmax)) &&
	g_xmax-g_xmin>0 &&  g_ymax-g_ymin >0 &&  gexp_zmax-gexp_zmin >0 ){

	  res.push_back(Box(g_xmin,g_ymin,gexp_zmin,g_xmax,g_ymax,gexp_zmax));
      }
      
    }

    void subtract(const Box& b1, const Box& b2, list<Box>& res, long64 pack_l, long64 pack_w){
      //b1 MUST have a non empty intersection w.r.t. the b2

      long64 b1_xmin= b1.getX();
      long64 b1_ymin= b1.getY();
      long64 b1_zmin= b1.getZ();
      long64 b2_xmin= b2.getX();
      long64 b2_ymin= b2.getY();
      long64 b2_zmin= b2.getZ();	
      long64 b1_xmax= b1.getX() + b1.getL();
      long64 b1_ymax= b1.getY() + b1.getW();
      long64 b1_zmax= b1.getZ() + b1.getH();
      long64 b2_xmax= b2.getX() + b2.getL();
      long64 b2_ymax= b2.getY() + b2.getW();
      long64 b2_zmax= b2.getZ() + b2.getH();
   
      if(b2_xmax<b1_xmax)
	    res.push_back(Box(b2_xmax,b1_ymin,b1_zmin,b1_xmax,b1_ymax,b1_zmax));
     
      if(b2_ymax<b1_ymax)
		res.push_back(Box(b1_xmin,b2_ymax,b1_zmin,b1_xmax,b1_ymax,b1_zmax));

      if(pack_l==-1 && b2_zmax<b1_zmax)
		res.push_back(Box(b1_xmin,b1_ymin,b2_zmax,b1_xmax,b1_ymax,b1_zmax));
	  else if(b2_zmax<b1_zmax)
		res.push_back(Box(b2_xmin,b2_ymin,b2_zmax,b2_xmin+pack_l,b2_ymin+pack_w,b1_zmax));
      
      if(b2_xmin>b1_xmin)
		res.push_back(Box(b1_xmin,b1_ymin,b1_zmin,b2_xmin,b1_ymax,b1_zmax));
      
      if(b2_ymin>b1_ymin)
		res.push_back(Box(b1_xmin,b1_ymin,b1_zmin,b1_xmax,b2_ymin,b1_zmax));
      
      if(b2_zmin>b1_zmin)
		res.push_back(Box(b1_xmin,b1_ymin,b1_zmin,b1_xmax,b1_ymax,b2_zmin));
      
    }

  
 
  bool Contain(btDbvtAabbMm& b, btDbvtAabbMm& a){
       return( (b.Mins().getX()<a.Mins().getX()+0.5)&&
                (b.Mins().getY()<a.Mins().getY()+0.5)&&
                 (b.Mins().getZ()<a.Mins().getZ()+0.5)&&
                 (b.Maxs().getX()+0.5>a.Maxs().getX())&&
                (b.Maxs().getY()+0.5>a.Maxs().getY())&&
                 (b.Maxs().getZ()+0.5>a.Maxs().getZ()));
  }

  void filter_nonmaximalGaps(list<Box>& gaps){
      //sorted gaps (e.g., y, luego z, luego x)
    list<Box>::iterator it=gaps.begin();

    for(;it!=gaps.end();it++){

	  list<Box>::iterator it2=it;
	  for(it2++;it2!=gaps.end();it2++){
	    if(Contain(*it,*it2)){
	      it2=gaps.erase(it2);
	      it2--;
	    }
	  }
    }
  }


  
  void print(const btDbvtAabbMm & b){
      printf("(%f,%f,%f);(%f,%f,%f)\n",b.Mins().getX(),b.Mins().getY(),b.Mins().getZ(),b.Maxs().getX(),b.Maxs().getY(),b.Maxs().getZ());
    }
  
  }
