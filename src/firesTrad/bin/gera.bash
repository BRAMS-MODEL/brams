#!/bin/bash
# ==============================================================================   
# title           : run_frp.sh
# description     : This script will get all files from GOES, MODIS, SEVERI, VIIRS and convert into FRE
# authors		      : GMAI Team
# date            : 20170201
# version         : 1.0    
# usage		 	  : bash create_emissions_frp_v.1.0.sh
# notes			  : Download folders goes-12, goes-13, goes-15
# ==============================================================================

# WORK DIR
SUBMIT_HOME=/dados/bramsrd/luiz.flavio/models/BRAMS-5.5/auxProgs/firesTrad/bin

#Local dos executaveis
BIN_HOME=/dados/bramsrd/luiz.flavio/models/BRAMS-5.5/auxProgs/firesTrad/bin

# LOCAL AND OUTPUT FOLDERS
dirProcessed=${SUBMIT_HOME}"processed/" ; mkdir -p ${SUBMIT_HOME}"processed/"
# Local onde estao os arquivos de entrada 
fire_data=/share/bramsrd/dist/BRAMS/data/EMISSIONS/
# Local onde estao os arquivos lulc
mdata=${SUBMIT_HOME}"data/"
# 
codes=${BIN_HOME}


# DATE IN FORMAT =  YYYY-MM-DD 
date_inic=2020-90-01


temp=${SUBMIT_HOME}"tmp/"   ; mkdir -p ${SUBMIT_HOME}"tmp/" 
# Numero de dias 
ndays=15

# Horas para processar nos arquivos da DSA
itime=0
ftime=24

# CONVERT FIRES TO 3BEM TRADITIONAL? YES or NO 
p3BEM="YES"

novo="NO"
#################### TO FRE ESTIMATION NOT FIRES TO 3BEM TRADITIONAL  #################################
# PROCESS FIRES YES or NO 
#DSA
pdsa='NO'

#MODIS EM HDF
pmodis='NO'

# MODIS TXT
tmodis='NO'

#GOES
pgoes='NO'

#VIIRS
pviirs='NO'

# From 2008 to 09:40 11-nov-2015 the METEOSAT is divided into regions            
# NAFR, SAFR, EURO,SAME to avoid to many files                                   
# From 09:45 11-nov-15 FULLDISK. 

#METEOSAT REGIONS
pmetereg='NO'

#METEOSAT FULL DISK
pmetefull='NO'

# LOCAL TO FEER DATA
feer=${mdata}"FEER.V1.bin"


#INPUT TO FIRES DATA_DIR 
#DSA - FILES IN DIR STRUCTURE YYYY/focos_dsa_YYYYMMDD.txt 
fires_dsa=${fire_data}"DSA/"

#MODIS MOD14 and MYD14 - FILES IN DIR STRUCTURE YYYY/DDD/MOD14.AYYYYDDD.HHMM.*.hdf 
fires_mod=${fire_data}"MOD14/"
fires_myd=${fire_data}"MYD14/"

#MODIS MOD14 and MYD14 ASCII
fires_tmod=${fire_data}"MODIS_V6_CCATT/TXT_GRANULE/MOD14/"
fires_tmyd=${fire_data}"MODIS_V6_CCATT/TXT_GRANULE/MYD14/"

#WFABBA GOES FIRE DATA - FILES IN DIR STRUCTURE YYYYDDD/fYYYYDDDHHMM.v65.g13.filt  
fires_wfa=${fire_data}"GOES16/"

#VIIRS FIRE DATA - FILES IN DIR STRUCTURE YYYYMMDD/AF_v1r1_npp_sYYYYMMDDHHMM586_e201711150013228_c201711150156380.txt
fires_vii=${fire_data}"VIIRS/"

# FULL YYYY/MM/HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_MSG-Disk_YYYYMMDDHHMM.bz2                                                                          
fires_met=${fire_data}"METEOSAT/" 

#AFRICA                                                                           
fires_metna=${fire_data}"METEOSAT/NAFRA/"
fires_metsa=${fire_data}"METEOSAT/SAFRA/"

#South America and Europe                                                                           
fires_meteu=${fire_data}"METEOSAT/EUROA/"
fires_metma=${fire_data}"METEOSAT/SAME/"

# LOCAL TO FRE CODES 
r_dsa=${codes}"dsa.exe"
r_goes=${codes}"goes.exe"
r_meteosat=${codes}"meteosat.exe"
r_modis=${codes}"modis.exe"
r_tmodis=${codes}"modis_txt.exe"
r_viirs=${codes}"viirs.exe"
r_frp=${codes}"frp.exe"
r_cfre=${codes}"fre.exe"

###################################################################
###################################################################
############## DO NOT CHANGE FROM HERE  ###########################
###################################################################
###################################################################

# Estimate the delta Day

delta_dayf=$(date -d "${date_inic}" +%j)
delta_day=$ndays

echo "********** INITIAL DATE:   ${date_inic} ******************"
echo "********** NUMBER OF DAYS: ${delta_day} ******************"

i=0
while [ ${i} -lt ${delta_day} ];
do
	data_proc=$(date -d  "${date_inic} ${i} days" +"%Y-%m-%d" )

echo "********** DATE: ${data_proc} ************************"
	
	dd=$(date -d  "$data_proc" +"%d" )
	mm=$(date -d  "$data_proc" +"%m" )
	yyyy=$(date -d  "$data_proc" +"%Y" )
	
	
 	juld=$(date -d "${yyyy}${mm}${dd}" +%j)
 	
 	rm -rf ${temp}/*

if [ ${novo} == "YES" ]; then


#############  DSA #################	

	if [ ${pdsa} == "YES" ]; then

		cd ${fires_dsa}${yyyy}${mm}

	  cp focos_dsa*${yyyy}${mm}${dd}*txt ${temp}

	  cd ${temp}
		
	      for file_dsa in focos_dsa*${yyyy}${mm}${dd}*txt
				     do
					       echo ${file_dsa}
				                  echo "PROCESSANDO ${r_dsa} ${file_dsa} ${itime} ${ftime} ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin ${yyyy}${juld}"  
					               				                 
	                        ${r_dsa} ${file_dsa} ${itime} ${ftime} ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin ${yyyy}${juld}  
	
													mv ${yyyy}${juld}* ${dirProcessed}
												
	      done
	fi

#############  FINISHED DSA #################		
	
#############  GOES  #################	

	if [ ${pgoes} == "YES" ]; then

		cd ${fires_wfa}${yyyy}/${juld}                                                                    
		                                                                                          
			cp f${yyyy}${juld}*.GOES-16 ${temp}

	    cd ${temp}
			
			for file_goes in f${yyyy}${juld}*.GOES-16                                                      
					do                                                                                  
							horag=$(echo ${file_goes} | cut -c9-12)                                             
			                                                                                        
		  				echo ${horag} ${file_goes}   
	  					                                                          
		          ${r_goes} ${file_goes} ${yyyy}${juld}${horag}goes.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horag} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin                 
		                                                                                          
						  mv ${yyyy}${juld}${horag}goes.out ${dirProcessed}                                   
		                                                                                          
		     done                                                                                 

	fi

#############  FINISHED GOES #################	

############# METEOSAT FULL DISK ##########################

	if [ ${pmetefull} == "YES" ]; then

		cd ${fires_met}${yyyy}/${mm}
	
	  cp *ListProduct*${yyyy}${mm}${dd}*bz2 ${temp}
	
	  cd ${temp}
	
			for file_met in HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_MSG-Disk_${yyyy}${mm}${dd}*.bz2
					do		
							bunzip2 ${file_met}
							fileHdf=$(echo ${file_met} | awk -F"." '{printf "%s", $1}')
							horam=$(echo ${file_met} | awk -F"_" '{printf "%s", $6}' | cut -c9-12 )
	
	            echo "PROCESSANDO ${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}metfull.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin"
	            
	            ${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}metfull.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin
	
						  mv ${yyyy}${juld}${horam}metfull.out ${dirProcessed}
	 
	      done  

	fi


########################## FINISHED METEOSAT FULL DISK ############################################

############# METEOSAT REGIONS ##########################
	
	if [ ${pmetereg} == "YES" ]; then

			cd ${fires_metna}${yyyy}/${mm}
		
		  cp *ListProduct*${yyyy}${mm}${dd}*bz2 ${temp}
		
		  cd ${temp}
		
				for filem in HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_NAfr_${yyyy}${mm}${dd}*.bz2
						do		
								bunzip2 ${filem}
								fileHdf=$(echo ${filem} | awk -F"." '{printf "%s", $1}')
								horam=$(echo ${filem} | awk -F"_" '{printf "%s", $6}' | cut -c9-12 )
		
								echo "PROCESSANDO ${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}metnaf.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin"
	            
	            	${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}metnaf.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin
		
							  mv ${yyyy}${juld}${horam}metnaf.out ${dirProcessed}
		 
		      done       
		
			cd ${fires_metsa}${yyyy}/${mm}
		
		  cp *ListProduct*${yyyy}${mm}${dd}*bz2 ${temp}
		
		  cd ${temp}
		
				for filem in HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_SAfr_${yyyy}${mm}${dd}*.bz2
						do		
								bunzip2 ${filem}
								fileHdf=$(echo ${filem} | awk -F"." '{printf "%s", $1}')
								horam=$(echo ${filem} | awk -F"_" '{printf "%s", $6}' | cut -c9-12 )
		
								echo "PROCESSANDO ${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}metsaf.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin"
	            
	            	${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}metsaf.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin
		
							  mv ${yyyy}${juld}${horam}metsaf.out ${dirProcessed}		
		      done  
		
			cd ${fires_meteu}${yyyy}/${mm}
		
		  cp *ListProduct*${yyyy}${mm}${dd}*bz2 ${temp}
		
		  cd ${temp}
		
				for filem in HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_Euro_${yyyy}${mm}${dd}*.bz2
						do		
								bunzip2 ${filem}
								fileHdf=$(echo ${filem} | awk -F"." '{printf "%s", $1}')
								horam=$(echo ${filem} | awk -F"_" '{printf "%s", $6}' | cut -c9-12 )
		
								echo "PROCESSANDO ${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}meteur.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin"
	            
	            	${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}meteur.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin
		
							  mv ${yyyy}${juld}${horam}meteur.out ${dirProcessed}
		 
		      done   
		
			cd ${fires_metma}${yyyy}/${mm}
		
		  cp *ListProduct*${yyyy}${mm}${dd}*bz2 ${temp}
		
		  cd ${temp}
		
				for filem in HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_SAme_${yyyy}${mm}${dd}*.bz2
						do		
								bunzip2 ${filem}
								fileHdf=$(echo ${filem} | awk -F"." '{printf "%s", $1}')
								horam=$(echo ${filem} | awk -F"_" '{printf "%s", $6}' | cut -c9-12 )
		
								echo "PROCESSANDO ${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}metsam.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin"
	            
	            	${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}metsam.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin
		
							  mv ${yyyy}${juld}${horam}metsam.out ${dirProcessed}
		 
		      done 

	fi
	
############# METEOSAT REGIONS FINISHED ##########################	
	
############# MODIS ##########################
	
	if [ ${pmodis} == "YES" ]; then	
	
			cd ${fires_mod}${yyyy}/${juld}                                                            
			 
			cp MOD*${yyyy}${juld}*hdf ${temp}

	    cd ${temp} 
			                                                                                          
			     for file_modis in MOD*${yyyy}${juld}*hdf                                                   
					     do                                                                               
						       echo ${file_modis}                                                                 
						            hora=$(echo ${file_modis} | awk -F"." '{printf "%s", $3}' )                   
			                       
			                       ${r_modis} ${file_modis} ${yyyy}${juld}${hora}mo.out ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${yyyy}${mm}${dd}${hora} ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin          
			                                                                                          
  														mv ${yyyy}${juld}${hora}mo.out ${dirProcessed}                          
			                                                                                          
			          done                                                                            
			                                                                                          
			cd ${fires_myd}${yyyy}/${juld}                                                            
			                                                                                          
			     cp MYD*${yyyy}${juld}*hdf ${temp} 

	         cd ${temp} 
			     
			     
			     for filemy in MYD*${yyyy}${juld}*hdf                                                 
					     do                                                                               
						       echo ${filemy}                                                               
						            horamy=$(echo ${filemy} | awk -F"." '{printf "%s", $3}' )               

			                       	${r_modis} ${filemy} ${yyyy}${juld}${horamy}my.out ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${yyyy}${mm}${dd}${horamy} ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin          
			                                                                                          
  														mv ${yyyy}${juld}${horamy}my.out ${dirProcessed}     
                  
			                                                                                          
			      		done	                                                                          
	
	fi	

############# MODIS FINISHED ##########################

############# VIIRS ##########################
	
	if [ ${pviirs} == "YES" ]; then	

			cd ${fires_vii}${yyyy}${mm}${dd} 
				
				cp AF_v1r2_npp_s${yyyy}${mm}${dd}*.txt ${temp}
				
			  cd ${temp}
			
					for filev in AF_v1r2_npp_s${yyyy}${mm}${dd}*.txt
							do
								
									horav=$(echo ${filev} | cut -c22-25)
					
					        echo ${horav} ${filev}
					        
			            ${r_viirs} ${filev} ${yyyy}${juld}${horav}viirs.out ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${yyyy}${mm}${dd}${horav} ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin

			 					  mv ${yyyy}${juld}${horav}viirs.out ${dirProcessed}
			 
			      done

	fi

############# MODIS TXT ##########################
	
	if [ ${tmodis} == "YES" ]; then	
	
			cd ${fires_tmod}${yyyy}/${juld}                                                            
			
			cp MOD*${yyyy}${juld}*txt ${temp}
				
			cd ${temp}
			                                                                                          
			     for file_tmodis in MOD*${yyyy}${juld}*txt                                               
					     do                                                                               
						       echo ${file_tmodis}                                                                 
						            horatm=$(echo ${file_tmodis} | awk -F"." '{printf "%s", $3}' )                   
			                       
			                       ${r_tmodis} ${file_tmodis} ${yyyy}${juld}${horatm}tmo.out ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${yyyy}${mm}${dd}${horatm} ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin          
			                                                                                          
  														mv ${yyyy}${juld}${horatm}tmo.out ${dirProcessed}                          
			                                                                                          
			          done                                                                            
			                                                                                          
			cd ${fires_tmyd}${yyyy}/${juld}                                                            
			                                                                                          
			cp MYD*${yyyy}${juld}*txt ${temp}
				
			cd ${temp}     
			     
			     
			     for file_tmy in MYD*${yyyy}${juld}*txt                                                
					     do                                                                               
						       echo ${file_tmy}                                                               
						            horatmy=$(echo ${file_tmy} | awk -F"." '{printf "%s", $3}' )               

			                       	${r_tmodis} ${file_tmy} ${yyyy}${juld}${horatmy}tmy.out ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${yyyy}${mm}${dd}${horatmy} ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin          
			                                                                                          
  														mv ${yyyy}${juld}${horatmy}tmy.out ${dirProcessed}     
                  
			                                                                                          
			      		done	                                                                          
	
	fi	

############# MODIS FINISHED ##########################

#############  INTEGRATION #################

    cd ${dirProcessed}                 
                                                                                
    cat ${yyyy}${juld}* > ${yyyy}${juld}.daily                                                 
    
    mkdir ${yyyy}${juld}
    
    mv  ${yyyy}${juld}.daily  ${yyyy}${juld} 
    
    rm -f *.out                                                                      
                                                                               
    cd ${dirProcessed}/${yyyy}${juld} 
    
    ${r_frp} ${yyyy}${juld}.daily ${yyyy}${juld}.frp            
    ${r_cfre} ${yyyy}${juld}.frp  f${yyyy}${juld}.oper.feer.filt ${feer} 
     

#    rm -f  ${yyyy}${juld}.daily ${yyyy}${juld}.frp  

     rm -rf ${temp}/*
       

############################################################################################################################################
fi

if [ ${p3BEM} == "YES" ]; then

dir3bem=${SUBMIT_HOME}"3BEM_processed/" ; mkdir -p ${SUBMIT_HOME}"3BEM_processed/"
d3bemdirgoes=${dir3bem}"GOES/" ; mkdir -p ${dir3bem}"GOES/"
d3bemdirdsa=${dir3bem}"DSA/" ; mkdir -p ${dir3bem}"DSA/"
d3bemdirMODIS=${dir3bem}"MODIS/" ; mkdir -p ${dir3bem}"MODIS/"
d3bemdirmeteosat=${dir3bem}"EXTRA/" ; mkdir -p ${dir3bem}"EXTRA/"

# LOCAL CODES 
d3bem_dsa=${codes}"dsa23bem.exe"
d3bem_goes=${codes}"goes23bem.exe"
d3bem_meteosat=${codes}"rmeteosat23bem.exe"
d3bem_modis=${codes}"rmodis23bem.exe"
d3bem_viirs=${codes}"rviirs23bem.exe"


#############  GOES  #################	
echo "Processando GOES *************************************"

		cd ${fires_wfa}${yyyy}/${juld}                                                                    
		cp f${yyyy}${juld}*.GOES-16 ${temp}

	  cd ${temp}
			
			for file_goes in f${yyyy}${juld}*.GOES-16                                                      
					do                                                                                  
							horag=$(echo ${file_goes} | cut -c9-12)                                             
			                                                                                        
		  				echo ${horag} ${file_goes}   
	  					                                                      
		          ${d3bem_goes} ${file_goes} f${yyyy}${juld}${horag}.v65.g16.filt ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horag} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin                 
		                                                                                          
						  mv f${yyyy}${juld}${horag}.v65.g16.filt ${d3bemdirgoes}                                   
		                                                                                          
		     done                                                                                 
      cd ${temp}
      rm -rf *


#############  FINISHED GOES #################	

############# METEOSAT FULL DISK ##########################
echo "Processando Meteosat ****************************"	

		cd ${fires_met}${yyyy}/${mm}
echo ${fires_met}${yyyy}/${mm}
	
	  cp *ListProduct*${yyyy}${mm}${dd}* ${temp}
	
	  cd ${temp}
ls
			for file_met in HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_MSG-Disk_${yyyy}${mm}${dd}*
					do		
							#bunzip2 ${file_met}
							fileHdf=$(echo ${file_met} | awk -F"." '{printf "%s", $1}')
							horam=$(echo ${file_met} | awk -F"_" '{printf "%s", $6}' | cut -c9-12 )
	
	            echo "PROCESSANDO ${r_meteosat} ${fileHdf} ${yyyy}${juld}${horam}metfull.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin"
	            
	            ${d3bem_meteosat} ${fileHdf} ${yyyy}${juld}${horam}metfull.out ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horam} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin
	
						  mv ${yyyy}${juld}${horam}metfull.out ${d3bemdirmeteosat}
           
	            
	      done  


							cd ${d3bemdirmeteosat}
	            
	            cat *metfull.out > Fires.extra.${yyyy}${mm}${dd}.txt
	            
	            rm -rf *.out

							cd ${temp}
							rm -rf *

########################## FINISHED METEOSAT FULL DISK ############################################

	
############# MODIS ##########################
#	
#			cd ${fires_mod}${yyyy}/${juld}                                                            
#			 
#			cp MOD*${yyyy}${juld}*hdf ${temp}
#
#	    cd ${temp} 
#			                                                                                          
#			     for file_modis in MOD*${yyyy}${juld}*hdf                                                   
#					     do                                                                               
#						       echo ${file_modis}                                                                 
#						            hora=$(echo ${file_modis} | awk -F"." '{printf "%s", $3}' )                   
#			                       
#			                       ${d3bem_modis} ${file_modis} ${yyyy}${juld}${hora}mo.out ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${yyyy}${mm}${dd}${hora} ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin          
#			                                                                                          
#  														mv ${yyyy}${juld}${hora}mo.out ${d3bemdirMODIS}                          
#			                                                                                          
#			          done                                                                            
#			                                                                                          
#			cd ${fires_myd}${yyyy}/${juld}                                                            
#			                                                                                          
#			     cp MYD*${yyyy}${juld}*hdf ${temp} 
#
#	         cd ${temp} 
#			     
#			     
#			     for filemy in MYD*${yyyy}${juld}*hdf                                                 
#					     do                                                                               
#						       echo ${filemy}                                                               
#						            horamy=$(echo ${filemy} | awk -F"." '{printf "%s", $3}' )               
#
#			                       	${d3bem_modis} ${filemy} ${yyyy}${juld}${horamy}my.out ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${yyyy}${mm}${dd}${horamy} ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin          
#			                                                                                          
 # 														mv ${yyyy}${juld}${horamy}my.out ${d3bemdirMODIS}      
  #                
#			                                                                                          
#			      		done	                                                                          
#	
#
############# MODIS FINISHED ##########################

############# VIIRS ##########################
echo "Processando VIIRS *********************************************"	
			cd ${fires_vii}${yyyy}${mm}${dd} 
				
				cp AF_v1r2_npp_s${yyyy}${mm}${dd}*.txt ${temp}
				
			  cd ${temp}
			
					for filev in AF_v1r2_npp_s${yyyy}${mm}${dd}*.txt
							do
								
									horav=$(echo ${filev} | cut -c22-25)
					
					        echo ${horav} ${filev}
					        
			            ${d3bem_viirs} ${filev} ${yyyy}${juld}${horav}viirs.out ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${yyyy}${mm}${dd}${horav} ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin

			 					  mv ${yyyy}${juld}${horav}viirs.out ${d3bemdirMODIS}
		 
			 
			      done


							cd ${d3bemdirMODIS}
	            
	            cat *.out > Fires.${yyyy}${mm}${dd}.txt
	            
	            rm -rf *.out
							cd ${temp}
							rm -rf *

#############  DSA #################	
echo "Processando DSA ***************************************"
		cd ${fires_dsa}

	  cp Focos*${yyyy}${mm}${dd}*txt ${temp}

	  cd ${temp}
		
	      for file_dsa in Focos*${yyyy}${mm}${dd}*txt
				     do
					       echo ${file_dsa}
				                  echo "PROCESSANDO ${r_dsa} ${file_dsa} ${itime} ${ftime} ${mdata}lulc/R2/${yyyy}.LUC.R2.bin ${mdata}R2/med/${yyyy}.MED.R2.bin ${mdata}R2/dvp/${yyyy}.DVP.R2.bin ${yyyy}${juld}"  
					               				                 
	                        ${d3bem_dsa} ${file_dsa} Focos.${yyyy}${mm}${dd}.proc.txt  
	
													mv Focos.${yyyy}${mm}${dd}.proc.txt ${d3bemdirdsa}
												
	      done
	      
	      cd ${temp}
				rm -rf *

#############  FINISHED DSA #################	





fi

((i=$i+1))
done





