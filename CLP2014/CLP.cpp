#include <iostream>
#include <fstream>
#include <sstream>
#include "CLP.h"
 
 using namespace odp;
 namespace odp{
	 CLP::CLP(ifstream& in, double total_time):total_time(total_time){
		 
		 begin_time=clock();
		 
		string line;
	 	getline(in, line ); //n_inst best_sol?	
	 	getline(in, line); //L W H
	 	std::stringstream ss(line);
	 	ss >> L >> W >> H;
	 	getline(in, line);
	 	std::stringstream ss0(line);
	 	int nb_types;
	 	ss0 >> nb_types;

	 	//se lee el archivo de entrada
	 	//Objetos BoxType, guardan los datos para cada tipo de cajas: dimensiones (w x l x h) 
	 	//y restricciones de rotación
	 	//En el objeto clp se agregan los tipos de cajas y el número de elementos que hay de cada tipo
	 	for(int j=0;j<nb_types;j++){
	 	  getline(in, line );
	 	  int n, id;
	 	  long64 l,h,w;
	 	  bool rot1, rot2, rot3;
	 	  std::stringstream ss1(line);
	 	  ss1 >> id >> l >> rot1 >> w >> rot2 >> h >> rot3 >> n;  
	 	  BoxType* boxt=new BoxType(l,w,h, id, rot1, rot2, rot3);
	 	  add_box_type(boxt,n);
	 	}
	 	
	 }
	 
	 //Genera un block para cada tipo de caja
	 set<Block*, compareBlocks> CLP::single_box_block_generator(){
	 	list<BoxType*>::iterator it=box_types.begin();
	 	set<Block*, compareBlocks> blocks;
	
	 	while(it!=box_types.end()){
	 		for(int o=0; o<6; o++){
	 			if((*it)->is_valid_orientation[(Orientation) o])
	 			     blocks.insert(new Block(*it,(Orientation) o)); 
	 		}
		
	 		it++;
	 	}
	 	return blocks;
	 }

	 //Ver paper
	 set<Block*, compareBlocks> CLP::general_block_generator(double min_fr, int max_bl,bool fsb){
	    clock_t begin=clock();	
	 	set<Block*, compareBlocks> B=single_box_block_generator();
		
		//set<Block*, compareBlocks>:: iterator it=B.begin();
		//for(;it!=B.end() ; it++)
			//(*it)->print();
		
	 	set<Block*, compareBlocks> P=B;
	 	while(B.size()<max_bl){
	          //  cout << B.size() << endl;
	 	    set<Block*, compareBlocks> N;
	 	    set<Block*, compareBlocks>:: iterator itP=P.begin();
		    int new_elems=0;
	 	    for(;itP!=P.end() && B.size()+new_elems<max_bl; itP++){
	 		set<Block*, compareBlocks>:: iterator itB=B.begin();
	 		for(;itB!=B.end() && B.size()+new_elems<max_bl ; itB++){	
					        
	 		   list<Block*> newB = Block::createBlocks(**itP, **itB, L,W,H,min_fr,fsb);
	 		   list<Block*>:: iterator itNew=newB.begin();
	 		   for(;itNew!=newB.end();itNew++){
	 		       if((*itNew)->feasible(nb_boxes)){
				     int NoldSize=N.size();
				     if(B.find(*itNew)==B.end() && N.find(*itNew)==N.end()){
				          N.insert(*itNew);   
				          new_elems++;			  
				     }else delete *itNew;
				     if(B.size()+new_elems>=max_bl) break;
			       }else delete *itNew;
	 		   }
	 	}
	  	        if(get_time()>20) return B; 
				
	 		}
	 		if(N.size()==0) break;
			B.insert(N.begin(),N.end());	
	 		P=N;
	 	}
	 	return B;
	 }
	 
   
 }
