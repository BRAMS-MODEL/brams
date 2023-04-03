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
SUBMIT_HOME=/dados/bramsrd/luiz.flavio/models/BRAMS-5.5/auxProgs/firesTrad/bin/

#Local dos executaveis
BIN_HOME=/dados/bramsrd/luiz.flavio/models/BRAMS-5.5/auxProgs/firesTrad/bin/
#/dados/bramsrd/luiz.flavio/models/BRAMS-5.5/auxProgs/firesTrad/bin/

# Local onde estao os arquivos de entrada 
fire_data=/share/bramsrd/dist/BRAMS/data/EMISSIONS/
out_data=/share/bramsrd/dist/BRAMS/data/FOCOS/
# Local onde estao os arquivos lulc
mdata=/dados/bramsrd/luiz.flavio/models/BRAMS-5.5/auxProgs/firesTrad/bin/data/

# 
codes=${BIN_HOME}


# DATE IN FORMAT =  YYYY-MM-DD 
date_inic=$1-$2-$3


temp=${SUBMIT_HOME}"tmp/"   ; mkdir -p ${SUBMIT_HOME}"tmp/" 
# Numero de dias 
ndays=1

# Horas para processar nos arquivos da DSA
itime=0
ftime=24

# CONVERT FIRES TO 3BEM TRADITIONAL? YES or NO 
p3BEM="YES"

# LOCAL TO FEER DATA
feer=${mdata}"FEER.V1.bin"


#INPUT TO FIRES DATA_DIR

#WFABBA GOES FIRE DATA - FILES IN DIR STRUCTURE YYYYDDD/fYYYYDDDHHMM.v65.g13.filt  
fires_wfa=${fire_data}"GOES16/"

#VIIRS FIRE DATA - FILES IN DIR STRUCTURE YYYYMMDD/AF_v1r1_npp_sYYYYMMDDHHMM586_e201711150013228_c201711150156380.txt
fires_vii=${fire_data}"VIIRS/"

# FULL YYYY/MM/HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_MSG-Disk_YYYYMMDDHHMM.bz2                                                                          
fires_met=${fire_data}"METEOSAT/" 

#DSA - FILES IN DIR STRUCTURE YYYY/focos_dsa_YYYYMMDD.txt 
fires_dsa=${fire_data}"DSA/"

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

dir3bem=${out_data}$1$2$3"/" ; mkdir ${dir3bem}
d3bemdirgoes=${dir3bem}"GOES/" ; mkdir -p ${dir3bem}"GOES/"
d3bemdirdsa=${dir3bem}"DSA/" ; mkdir -p ${dir3bem}"DSA/"
d3bemdirMODIS=${dir3bem}"MODIS/" ; mkdir -p ${dir3bem}"MODIS/"
d3bemdirmeteosat=${dir3bem}"EXTRA/" ; mkdir -p ${dir3bem}"EXTRA/"

echo "::::: 3bem  dir ::::: ${dir3bem} :::::: "
echo "::::: GOES  dir ::::: ${d3bemdirgoes} :::::: "
echo "::::: DSA   dir ::::: ${d3bemdirdsa} :::::: "
echo "::::: MODIS dir ::::: ${d3bemdirMODIS} :::::: "
echo "::::: METEOSAT dir ::::: ${d3bemdirmeteosat} :::::: "

	
# LOCAL CODES 
d3bem_dsa=${codes}"dsa23bem.exe"
d3bem_goes=${codes}"goes23bem.exe"
d3bem_meteosat=${codes}"rmeteosat23bem.exe"
d3bem_modis=${codes}"rmodis23bem.exe"
d3bem_viirs=${codes}"rviirs23bem.exe"

i=0
while [ ${i} -lt ${delta_day} ];
do
	data_proc=$(date -d  "${date_inic} ${i} days" +"%Y-%m-%d" )
	echo "********** DATE: ${data_proc} ************************"
#
	dd=$(date -d  "$data_proc" +"%d" )
	mm=$(date -d  "$data_proc" +"%m" )
	yyyy=$(date -d  "$data_proc" +"%Y" )
	juld=$(date -d "${yyyy}${mm}${dd}" +%j)
	rm -rf ${temp}/*
#	
	#############  GOES  #################	
	if [ $4 -gt 0 ]
	then
	echo "Processando GOES *************************************"
		cd ${fires_wfa}${yyyy}/${juld}                                                                    
		cp f${yyyy}${juld}*.FDCF.GOES-16 ${temp}
	
		cd ${temp}
		
		for file_goes in f${yyyy}${juld}*.GOES-16                                                      
		do                                                                                  
			horag=$(echo ${file_goes} | cut -c9-12)                                             
			hora2=$(echo ${file_goes} | cut -c11-12)
			if [ $hora2 -eq "00" ] || [  $hora2 -eq "30" ]
			then
				echo ${horag} ${file_goes}   
				${d3bem_goes} ${file_goes} f${yyyy}${juld}${horag}.v65.g16.filt ${mdata}lulc/R7/${yyyy}.LUC.R7.bin ${yyyy}${mm}${dd}${horag} ${mdata}R7/med/${yyyy}.MED.R7.bin ${mdata}R7/dvp/${yyyy}.DVP.R7.bin                 
				mv f${yyyy}${juld}${horag}.v65.g16.filt ${d3bemdirgoes}                                   
			fi
		done                                                                                 
	    cd ${temp}
	    rm -rf *
	#############  FINISHED GOES #################	
	fi
	
	############# METEOSAT FULL DISK ##########################
	if [ $5 -gt 0 ]
	then
		echo "Processando Meteosat ****************************"	
		
		cd ${fires_met}${yyyy}/${mm}/
		echo ${fires_met}${yyyy}/${mm}/
		cp *ListProduct*${yyyy}${mm}${dd}* ${temp}
		cd ${temp}
	        bunzip2 *.bz2
		for file_met in HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_MSG-Disk_${yyyy}${mm}${dd}*
		do		
			bunzip2 ${file_met}
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
	fi
	
	############# VIIRS ##########################
	if [ $6 -gt 0 ]
	then
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
	fi
	
	#############  DSA #################	
	if [ $7 -gt 0 ]
	then
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
