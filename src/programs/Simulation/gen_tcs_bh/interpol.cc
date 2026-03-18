#include "interpol.h"
 using namespace std;  

//-------------------------------------------------------------------------------------
double linear_interpol_tcs6D_4(
        TableProvider ddvcs_table_provider,
        float &cross_BH, float &cross_TCS, float &polpp,
        float &polpm , float &polmp , float &polmm , float &bhpara, float &bhperp,
        double mt, double Q2, double PhiCM, double ThetaCM,
        float mt_grid[],float Q2_grid[], float ThetaCM_grid[], float PhiCM_grid[],
	int targetpoldir, int beampoltype,
        int i, int j, int k, int n, int o,int l
){
	const bool debug=0;
        double result[8]={0};
        double f[8][2][2][2][2] ;
        Vars_Sef_tcs_pol sef[2][2][2][2] ;
	int jb, kb, nb, ob;

	if (debug){
                std::cout<<" "<<std::endl;
                std::cout << "[DT] Table size after FillTable() : " << ddvcs_table_provider.size() << std::endl;
		std::cout<<"bins: "<<i<<" " <<j<<" "<<k<<" "<<n<<" "<<o<<" "<<l<<std::endl;

                std::cout<<"var= "<<" "<<mt<<" "<<Q2<< " "<<PhiCM<<" "<<ThetaCM<<" "<<std::endl;
                std::cout<<"lm= "<<" "<<mt_grid[j]<<" "<<Q2_grid[k]<<" "<<PhiCM_grid[o]*180./3.1415<<" "<<ThetaCM_grid[n]*180./3.1415<<" "<<std::endl;
                std::cout<<"lM= "<<" "<<mt_grid[j+1]<<" "<<Q2_grid[k+1]<<" "<<PhiCM_grid[o+1]*180./3.1415<<" "<<ThetaCM_grid[n+1]*180./3.1415<<std::endl;
                std::cout<<"check angles (phi, th): "<<PhiCM*180./3.1415<<" "<<ThetaCM*180./3.1415<<std::endl;
                std::cout<<"example from the table: "<< ddvcs_table_provider.get_sect_eff_tcs_pol((Vars_Cin_tcs_pol){{i,j,k,n,o,l}}) <<std::endl;
                Vars_Sef_tcs_pol fs= ddvcs_table_provider.get_sect_eff_tcs_pol((Vars_Cin_tcs_pol){{i,j,k,n,o,l}});
                for (int p=0;p<6;p++){
                        cout<<"entry "<<p<<" : "<<getval_tcs_pol( fs  ,p)<< " ";
                } cout<<" "<<endl;
        }

                for (int jj=0;jj<2;jj++){
                        for (int kk=0;kk<2;kk++){
                                for (int nn=0;nn<2;nn++){
                                        for (int oo=0;oo<2;oo++){
						if (j+jj<NT) jb=j+jj; else jb=j; 
						if (k+kk<NQp2) kb=k+kk; else kb=k; 
						if (n+nn<NTh) nb=n+nn; else nb=n; 
						if (o+oo<NPhi) ob=o+oo; else ob=0; 

						sef[jj][kk][nn][oo] = ddvcs_table_provider.get_sect_eff_tcs_pol((Vars_Cin_tcs_pol){{i, jb, kb, nb, ob, l}});
						for  (int crA=0;crA<8; crA++){
							if (beampoltype==0 && crA>=6) break;
							f[crA][jj][kk][nn][oo]=getval_tcs_pol(sef[jj][kk][nn][oo],crA); 
						}

	}}}}


	for (int crA=0;crA<8; crA++){

		if (beampoltype==0 && crA>=6) break;
		if (beampoltype==1 && crA==0) continue;
                   for (int nn=0;nn<2;nn++){
                     for (int oo=0;oo<2;oo++){
            		    for (int jj=0;jj<2;jj++){
                        	for (int kk=0;kk<2;kk++){

				if (f[crA][jj][kk][nn][oo]==0){
					if (debug==1) cout<<"\n ref bin at 0 "<<" "<<jj<<" "<<kk<<" ";	
					if (jj==0){ 
						f[crA][jj][kk][nn][oo]=f[crA][1][0][nn][oo];
					} else if (jj==1 && kk==1){ 
						if (f[crA][0][0][nn][oo]!=0) f[crA][jj][kk][nn][oo]=f[crA][0][0][nn][oo];
					} else { 
						if (f[crA][1][0][nn][oo]!=0) f[crA][jj][kk][nn][oo]=f[crA][1][0][nn][oo];
						else f[crA][jj][kk][nn][oo]= f[crA][1][1][nn][oo]; 
					}
					if (debug==1) cout<<f[crA][jj][kk][nn][oo]<<endl;
				}
	    }}}}


	    result[crA] =	
		1./((mt_grid[j+1]-mt_grid[j] ) * (Q2_grid[k+1]-Q2_grid[k] ) * (PhiCM_grid[o+1]-PhiCM_grid[o] ) * (ThetaCM_grid[n+1]-ThetaCM_grid[n] ) )
	   *(	
		   (mt_grid[j+1]-mt) *(
		      (Q2_grid[k+1]-Q2) *(
	        	 (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) *
					  f[crA][0][0][0][0]
			    	+(ThetaCM-ThetaCM_grid[n]) *
					f[crA][0][0][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) *
                                      f[crA][0][0][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) *
                                      f[crA][0][0][1][1]
	      		   )
      	    		  )
	        	+ (Q2-Q2_grid[k] ) * (
	        	   (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) *
                                        f[crA][0][1][0][0]
			    +(ThetaCM-ThetaCM_grid[n])*
                                        f[crA][0][1][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) *
                                       f[crA][0][1][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) *
                                       f[crA][0][1][1][1]
	      		   )
	        	)
		   )
		   + (mt-mt_grid[j] ) *(
		      (Q2_grid[k+1]-Q2) *(
	        	 (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) *
                                        f[crA][1][0][0][0]
			    +(ThetaCM-ThetaCM_grid[n])*
                                        f[crA][1][0][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM)*
                                       f[crA][1][0][1][0]
			     +(ThetaCM-ThetaCM_grid[n])*
                                       f[crA][1][0][1][1]
	      		   )
      	    		  )
	        	+ (Q2-Q2_grid[k] ) * (
	        	   (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) *
                                       f[crA][1][1][0][0]
			    +(ThetaCM-ThetaCM_grid[n])*
                                       f[crA][1][1][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) *
                                       f[crA][1][1][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) *
                                       f[crA][1][1][1][1]
	      		   )
	        	)
		)
            );
	}

	if (beampoltype==0) cross_BH = result[0];
	else cross_BH = 0.5*(result[6]+result[7]);
	cross_TCS = result[1];
	polpp = result[2];
	polpm = result[3];
	polmp = result[4];
	polmm = result[5];
	if (beampoltype==1){
		bhpara=result[6]; 
		bhperp=result[7]; 
	} else {
		bhpara=0; 
		bhperp=0;
	}
	if (debug) {
		cout<<"res "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<result[4]<<" "<<result[0]<<endl;
		cout<<"pol "<<polpp<<" " <<polpm<<" "<<polmp<<" "<<polmm<<endl;
		cout<<"diff "<<cross_BH<<" "<<cross_TCS<<" "<<bhpara<<" "<<bhperp<<endl;
	}	
	return 0.25*(polpp+polmp+polmm+polpm);

} // end 6D function interpol 4


//-------------------------------------------------------------------------------------

double linear_interpol_tcs5D( TableProvider ddvcs_table_provider,float &cross_BH,float &cross_TCS, 
	float &polpp, float &polpm, float &polmp, float &polmm,
	double xbj, double mt, double Q2, double PhiCM, double ThetaCM,
	float xbj_grid[], float mt_grid[],float Q2_grid[], float ThetaCM_grid[], float PhiCM_grid[],
	int i, int j, int k, int n, int o){

	const bool debug=0; 
	double result[6]={0};
	double f[6][2][2][2][2][2] ;
	Vars_Sef_tcs sef[2][2][2][2][2] ; 
	int ib, jb, kb, nb, ob;

	if (debug){
		std::cout<<" "<<std::endl;
		std::cout << "[DT] Table size after FillTable() : " << ddvcs_table_provider.size() << std::endl;  
		std::cout<<"var= "<<xbj<<" "<<mt<<" "<<Q2<< " "<<PhiCM<<" "<<ThetaCM<<std::endl;  
		std::cout<<"lm= "<<xbj_grid[i]<<" "<<mt_grid[j]<<" "<<Q2_grid[k]<<" "<<PhiCM_grid[o]*180./3.1415<<" "<<ThetaCM_grid[n]*180./3.1415<<std::endl;  
		std::cout<<"lM= "<<xbj_grid[i+1]<<" "<<mt_grid[j+1]<<" "<<Q2_grid[k+1]<<" "<<PhiCM_grid[o+1]*180./3.1415<<" "<<ThetaCM_grid[n+1]*180./3.1415<<std::endl; 
		std::cout<<"check angles (phi, th): "<<PhiCM*180./3.1415<<" "<<ThetaCM*180./3.1415<<std::endl; 
		std::cout<<"example from the table: "<< ddvcs_table_provider.get_sect_eff_tcs((Vars_Cin_tcs){{i,j,k,n,o}}) <<std::endl;
		Vars_Sef_tcs fs= ddvcs_table_provider.get_sect_eff_tcs((Vars_Cin_tcs){{i,j,k,n,o}});
		for (int p=0;p<6;p++){
			cout<<"entry "<<p<<" : "<<getval_tcs( fs  ,p)<< " ";
		}cout<<" "<<endl;
	}

	  for (int ii=0;ii<2;ii++){
	  	for (int jj=0;jj<2;jj++){
	  		for (int kk=0;kk<2;kk++){
				for (int nn=0;nn<2;nn++){
					for (int oo=0;oo<2;oo++){
						if (i+ii<NEB) ib=i+ii; else ib=i;	
						if (j+jj<NT) jb=j+jj; else jb=j; 
						if (k+kk<NQp2) kb=k+kk; else kb=k; 
						if (n+nn<NTh) nb=n+nn; else nb=n; 
						if (o+oo<NPhi) ob=o+oo; else ob=0; 

						sef[ii][jj][kk][nn][oo] = ddvcs_table_provider.get_sect_eff_tcs((Vars_Cin_tcs){{ib, jb, kb, nb, ob}});
						if (debug==1){
							cout<<ii<<" "<<jj<<" "<<kk<<" "<<nn<<" "<<oo<<" / "<<ib<<" "<<jb<<" "<<kb<<" "<<nb<<" "<<ob<<" "<< sef[ii][jj][kk][nn][oo]  << endl;
						}

						for (int crA=0;crA<6; crA++){
							f[crA][ii][jj][kk][nn][oo]=getval_tcs(sef[ii][jj][kk][nn][oo],crA);
						}


	  }}}}}

    for (int crA=0;crA<6; crA++){
	for (int ii=0;ii<2;ii++){
	   for (int nn=0;nn<2;nn++){
		for (int oo=0;oo<2;oo++){
	  	    for (int jj=0;jj<2;jj++){
	  		for (int kk=0;kk<2;kk++){

				if (f[crA][ii][jj][kk][nn][oo]==0){
					if (debug==1) cout<<"\n fffffffffffff0 "<<ii<<" "<<jj<<" "<<kk<<" ";	
					if (jj==0){ // tmin
						if (f[crA][ii][1][0][nn][oo]!=0) f[crA][ii][jj][kk][nn][oo]=f[crA][ii][1][0][nn][oo];
						else if (f[crA][1][1][0][nn][oo]!=0)	f[crA][ii][jj][kk][nn][oo]=f[crA][1][1][0][nn][oo];
					} else if (jj==1 && kk==1){ // Q2max
						if (f[crA][ii][0][0][nn][oo]!=0) f[crA][ii][jj][kk][nn][oo]=f[crA][ii][0][0][nn][oo];
						else if (f[crA][1][0][0][nn][oo]!=0) f[crA][ii][jj][kk][nn][oo]=f[crA][1][0][0][nn][oo];
					} else { // jj=1 kk=0
						if (f[crA][1][1][0][nn][oo]!=0) f[crA][ii][jj][kk][nn][oo]=f[crA][1][1][0][nn][oo];
						else f[crA][ii][jj][kk][nn][oo]= f[crA][1][1][1][nn][oo]; 
					}
					if (debug==1) cout<<f[crA][ii][jj][kk][nn][oo]<<endl;
				}

	}}}}}


	   result[crA] = 
	  1./((xbj_grid[i+1]-xbj_grid[i] ) *(mt_grid[j+1]-mt_grid[j] ) * (Q2_grid[k+1]-Q2_grid[k] ) * (PhiCM_grid[o+1]-PhiCM_grid[o] ) * (ThetaCM_grid[n+1]-ThetaCM_grid[n] ) ) 
	  *(
	       (xbj_grid[i+1]-xbj) * (
		   (mt_grid[j+1]-mt) *(
		      (Q2_grid[k+1]-Q2) *(
	        	 (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) * f[crA][0][0][0][0][0]
			    	+(ThetaCM-ThetaCM_grid[n])* f[crA][0][0][0][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) * f[crA][0][0][0][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) * f[crA][0][0][0][1][1]
	      		   )
      	    		  )
	        	+ (Q2-Q2_grid[k] ) * (
	        	   (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) * f[crA][0][0][1][0][0]
			    +(ThetaCM-ThetaCM_grid[n]) * f[crA][0][0][1][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) * f[crA][0][0][1][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) * f[crA][0][0][1][1][1]
	      		   )
	        	)
		   )
		   + (mt-mt_grid[j] ) *(
		      (Q2_grid[k+1]-Q2) *(
	        	 (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) * f[crA][0][1][0][0][0]
			    +(ThetaCM-ThetaCM_grid[n]) *    f[crA][0][1][0][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) * f[crA][0][1][0][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) * f[crA][0][1][0][1][1]
	      		   )
      	    		  )
	        	+ (Q2-Q2_grid[k] ) * (
	        	   (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) * f[crA][0][1][1][0][0]
			    +(ThetaCM-ThetaCM_grid[n]) *    f[crA][0][1][1][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) * f[crA][0][1][1][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) * f[crA][0][1][1][1][1]
	      		   )
	        	)
		   )
        	  )
	       + (xbj-xbj_grid[i] )* (
		   (mt_grid[j+1]-mt) *(
		      (Q2_grid[k+1]-Q2) *(
	        	 (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) * f[crA][1][0][0][0][0]
			    +(ThetaCM-ThetaCM_grid[n])    * f[crA][1][0][0][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) * f[crA][1][0][0][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) * f[crA][1][0][0][1][1]
	      		   )
      	    		  )
	        	+ (Q2-Q2_grid[k] ) * (
	        	   (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) * f[crA][1][0][1][0][0]
			    +(ThetaCM-ThetaCM_grid[n]) *    f[crA][1][0][1][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) * f[crA][1][0][1][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) * f[crA][1][0][1][1][1]
	      		   )
	        	)
		   )
		   + (mt-mt_grid[j] ) *(
		      (Q2_grid[k+1]-Q2) *(
	        	 (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) * f[crA][1][1][0][0][0]
			    +(ThetaCM-ThetaCM_grid[n]) *    f[crA][1][1][0][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) * f[crA][1][1][0][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) * f[crA][1][1][0][1][1]
	      		   )
      	    		  )
	        	+ (Q2-Q2_grid[k] ) * (
	        	   (PhiCM_grid[o+1]-PhiCM) *(
	      	        	(ThetaCM_grid[n+1]-ThetaCM) * f[crA][1][1][1][0][0]
			    +(ThetaCM-ThetaCM_grid[n]) *    f[crA][1][1][1][0][1]
	      		   )
			   + (PhiCM-PhiCM_grid[o])*(
	      		      (ThetaCM_grid[n+1]-ThetaCM) * f[crA][1][1][1][1][0]
			     +(ThetaCM-ThetaCM_grid[n]) * f[crA][1][1][1][1][1]
	      		   )
	        	)
		   )
        	  )
	    );

	}

	cross_TCS=result[1];
	cross_BH=result[0];
	polpp=result[2];
        polpm = result[3];
        polmp = result[4];
        polmm = result[5];

	return (result[2]+result[3]+result[4]+result[5])*0.25;

}

// end interpolate

///////////////////////////////////////::

double getval_tcs(Vars_Sef_tcs sef,int crA){

	if (crA==0){
		return sef.vars[Vars_Sef_tcs::BH];
	} else if (crA==1){
		return sef.vars[Vars_Sef_tcs::TCS];
	} else if (crA==2){
		return sef.vars[Vars_Sef_tcs::SUMPP];
	} else if (crA==3){
		return sef.vars[Vars_Sef_tcs::SUMPM];
	} else if (crA==4){
		return sef.vars[Vars_Sef_tcs::SUMMP];
	} else if (crA==5){
		return sef.vars[Vars_Sef_tcs::SUMMM];
	} else return 0;

}

///////////////////////////////////////////////
double getval_tcs_pol(Vars_Sef_tcs_pol sef,int crA){

	if (crA==0){
		return sef.vars[Vars_Sef_tcs_pol::BH];
	} else if (crA==1){
		return sef.vars[Vars_Sef_tcs_pol::TCS];
	} else if (crA==2){
		return sef.vars[Vars_Sef_tcs_pol::SUMPP];
	} else if (crA==3){
		return sef.vars[Vars_Sef_tcs_pol::SUMPM];
	} else if (crA==4){
		return sef.vars[Vars_Sef_tcs_pol::SUMMP];
	} else if (crA==5){
		return sef.vars[Vars_Sef_tcs_pol::SUMMM];
	} else if (crA==6){
		return sef.vars[Vars_Sef_tcs_pol::BHpara];
	} else if (crA==7){
		return sef.vars[Vars_Sef_tcs_pol::BHperp];
	} else return 0;

}




