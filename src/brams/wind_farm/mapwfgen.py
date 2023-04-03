#!/usr/bin/env python
"""
@package BRAMS-WindFram
@author Luiz Flavio
@date 06/07/2017
@brief Script to read and plot a map from RAMSIN & Turbine files of BRAMS windfarm module
@copyright Under CC-GPL license
@link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
@param see the help (-h)
@version 1.0
@warning Make sure that numpy, h5py, matplotlib,python-tk,pygrads,pylab and windrose are installed
@warning sudo pip install h5py;sudo pip install mapplotlib;sudo pip install python-tk;sudo pip install windrose
"""
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import sys, getopt
import grads
from grads.ganum import GaNum 
import numpy as np
from numpy import meshgrid
from pylab import contourf
from numpy import linspace
import h5py
import datetime
import calendar
import csv
from windrose import WindroseAxes

def main(argv):

	err=printHeader()

	try:
		opts, args = getopt.getopt(argv,"hn:w:g:i:a:c:r:R:t:f:T:o:S",
			["ramsin=","wtable=","gradsf=","ihdf5=","arrow=","city=","roads=","rivers=",
			"turbines=","fillcont=","title=","output=","show="])
	except getopt.GetoptError:
		err=printHelp()
		sys.exit(2)

	#Reading the arguments from system call
	(ramsin,table,gradsf,ihdf5,plotarrow,plotcity,plotroads,
	plotrivers,plotturbines,fillcont,ptitle,output,show)=readArgs(opts, args)


	#readind the RAMSIN
	(totx,toty,centlat,centlon)=ReadRamsin(ramsin)
	
	if table!='':
		#reading table of turbines
		(xt,yt,windfarm,turb,nturb)=readTurb(table)
	
	#Making the basic map
	m=doBasicMap(totx,toty,centlat,centlon)
	
	if gradsf!='':
		#Opening and reading grads file 
		(ga,xa,ya)=readGrads(gradsf,m)
	
	fig1=plt.figure(1)
	fig1.set_size_inches(18.5, 10.5, forward=True)
	plt.suptitle(ptitle)


	plt.subplot(1,2,1)
	if gradsf!='':
		#Ploting Temp and Topo at window 1
		err=plotTempkTopo(ga,xa,ya,m)	
	#Putting something more on maps
	err=cherry(m,plotcity,plotroads,plotrivers,fillcont)
	
	plt.subplot(1,2,2)
	if gradsf!='':
		#Ploting winds at window 2
		err=plotWinds(ga,xa,ya,m,plotarrow)
	#Putting something more on maps
	err=cherry(m,plotcity,plotroads,plotrivers,fillcont)
	
	if table!='':
		#Putting marks inside map with turbines position
		err=plotTurb(plotturbines,xt,yt,m,turb,nturb)

	fig1.savefig(output+'.png',dpi=300)

	if ihdf5!='':
		(power,tini,tend,u,v,wind,dir,t,windfarm,turbine,pgeo,hdf5Info)=readHdf5(ihdf5)
		
		#Writing a csv file with turbine data 
		err=writeCSV(ihdf5,power,tini,u,v,wind,dir,hdf5Info)
	
		#####################################
		#Generating a graphic with power & winds
		#####################################
		fig2=plt.figure(2)
		fig2.set_size_inches(18.5, 10.5, forward=True)
		plt.suptitle(ptitle)
	
		#Plotting the power 
		plt.subplot(211)
		err=plotPower(t,power,hdf5Info,windfarm,turbine,pgeo,tini,tend)
	
		#plotting the winds
		plt.subplot(212)
		err=plotuv(t,u,v,wind,hdf5Info)
		fig2.savefig(windfarm[0]+'#'+turbine[0]+'.png',dpi=300)
		
		#####################################
		#Generating a graphic with windrose
		#####################################	
		#plotting the windrose
		plotWindRose(dir,wind,windfarm,turbine,pgeo,tini,tend)

		#writing a figure file with windrose
		plt.savefig(windfarm[0]+'#'+turbine[0]+'-WR.png',dpi=300)
	
	if(show=='1'):
		plt.show()




###############################################################################
################################ functions ####################################

def printHelp():
	print ''
	print 'geraFigs.py -n <ramsinFile> -w <windfarmsFile> -g <gradsFile (*1) -i <turbine_hdf5 (*1)> --output <output (*1)>'
	print 'Notes: '
	print '  # (*1) with no extension (ctl, bin, png or h5)'
	print '  # Other options:'
	print '             -a <0|1> 1 for wind vectors, 0 for streamlines  (--arrow)'
	print '             -c <0|1> 1 plots towns and cities of Brazil (--city)' 
	print '             -r <0|1> 1 plots roads of Brazil (--roads)'
	print '             -R <0|1> 1 plots rivers of Brazil (--rivers)'
	print '             -t <0|1> 1 plots the turbines position (--turbines)'
	print '             -f <0|1> 1 fills the continental area (*2) (--fillcont)'
	print '             -T Title for the plotting (--title)'
	print '             --output Name of output figure'
	print '             -S <0|1> 0 no show on screen, only write files (--show)'
	print '  # (*2) The option -f hides the others plots'
	print ''
	print 'Example: '
	print 'geraFigs.py -n RAMSIN-5km -w wturbines.txt -g os-output-g1 -a 0 -t 1 -T "WindFarm results" --output "wt"'
	return

def printHeader():
	print ''
	print ''
	print "                           ''~``"
	print "                          ( o o )"
	print "+--------------------.oooO--(_)--Oooo.-------------------+"
	print "|         Python Script to plot map & vars from BRAMS    |"
	print "| BRAMS - Brasilian Regional Atmospheric Modeling System |"
	print "|                                                        |"
	print "|    by: Luiz Flavio Rodrigues INPE/CPTEC   May 2017     |"
	print "|                                                        |"
	print "|                      .oooO                             |"
	print "|                      (   )   Oooo.                     |"
	print "+-----------------------\ (----(   )---------------------+"
	print "                         \_)    ) /"
	print "                               (_/"
	print ""
	print ""

def readArgs(opts, args):
	ramsin='RAMSIN'      #-n
	table=''             #-w
	gradsf=''            #-g
	ihdf5=''             #-i
	plotarrow='1'        #-a
	plotcity='1'         #-c
	plotroads='0'        #-r
	plotrivers='0'       #-R
	plotturbines='0'     #-t
	fillcont='0'         #-f
	ptitle=''            #-T
	output='saida'       #-F
	show='0'             #-S
	for opt, arg in opts:
		if opt == '-h':
			err=printHelp()
			sys.exit()
		elif opt in ("-n", "--ramsin"):
			ramsin = arg
		elif opt in ("-w", "--wtable"):
			table = arg
		elif opt in ("-g", "--gradsf"):
			gradsf = arg
		elif opt in ("-i", "--ihdf5"):
			ihdf5 = arg
		elif opt in ("-a", "--arrow"):
			if arg!='0' and arg!='1':
				print arg
				print 'argument for "-a" must be 0 or 1'
				sys.exit(2)
			plotarrow = arg
		elif opt in ("-c", "--city"):
			if arg!='0' and arg!='1':
				print 'argument for "-c" must be 0 or 1'
				sys.exit(2)
			plotcity = arg
		elif opt in ("-r", "--roads"):
			if arg!='0' and arg!='1':
				print 'argument for "-r" must be 0 or 1'
				sys.exit(2)
			plotroads = arg
		elif opt in ("-R", "--rivers"):
			if arg!='0' and arg!='1':
				print 'argument for "-R" must be 0 or 1'
				sys.exit(2)
			plotrivers = arg
		elif opt in ("-t", "--turbines"):
			if arg!='0' and arg!='1':
				print 'argument for "-t" must be 0 or 1'
				sys.exit(2)
			plotturbines = arg
		elif opt in ("-f", "--fillcont"):
			if arg!='0' and arg!='1':
				print 'argument for "-f" must be 0 or 1'
				sys.exit(2)
			fillcont = arg
		elif opt in ("-T", "--title"):
			ptitle = arg
		elif opt in ("-o", "--output"):
			output = arg	
		elif opt in ("-S", "--show"):
			if arg!='0' and arg!='1':
				print 'argument for "-S" must be 0 or 1'
				sys.exit(2)
			show=arg

	return ramsin,table,gradsf,ihdf5,plotarrow,plotcity,plotroads,plotrivers,plotturbines,fillcont,ptitle,output,show

def ReadRamsin(ramsin):
	#Try open ranmsin file send by argument #1
	try:
		f=open(ramsin,'r')
	except:
		print 'ERROR: RAMSIN not found. Check if file '+ramsin+' is available'
		sys.exit(0)
	
	line=f.readlines()
	f.close

	i=0
	posi=0
	while i < len(line):
		data=line[i].split()
		if len(data)>2:
			if data[0]=='NNXP':
				nnxp=float(data[2].replace(",", ""))
			if data[0]=='NNYP':
				nnyp=float(data[2].replace(",", ""))
			if data[0]=='DELTAX':
				deltax=float(data[2].replace(",", ""))
			if data[0]=='DELTAY':
				deltay=float(data[2].replace(",", ""))
			if data[0]=='CENTLAT':
				centlat=float(data[2].replace(",", ""))
			if data[0]=='CENTLON':
				centlon=float(data[2].replace(",", ""))
		i=i+1
	#Calculating the distance in each direction
	totx=nnxp*deltax
	toty=nnyp*deltay
	return totx,toty,centlat,centlon

def readTurb(table):
	xt=[]
	yt=[]
	turb=[]
	try:
		f=open(table,'r')
	except:
		print 'ERROR: Wind table not found. Check if file '+table+' is available'
		sys.exit(0)
	
	line=f.readlines()
	f.close
	i=0
	nturb=len(line)
	#getting turbine position and windfarm/turbine name
	while i < len(line):
		data=line[i].split()
		xt.append(float(data[1]))
		yt.append(float(data[0]))
		windfarm=data[3]
		turb.append(data[4])
		i=i+1

	return xt,yt,windfarm,turb,nturb

def cherry(m,plotcity,plotroads,plotrivers,fillcont):
	m.drawcountries()
	m.drawstates()
	#Drawing cities and towns 
	if plotcity=='1':
		m.readshapefile('/home/lufla/workspace/brams-wind/Shapefiles/municipios_2010', 'cidades')
	#Drawing roads
	if plotroads=='1': 
		m.readshapefile('/home/lufla/workspace/brams-wind/Shapefiles/ST_DNIT_Rodovias_SNV2015_03', 'rodovias',color='m')
	#Drawing rivers
	if plotrivers=='1':
		m.readshapefile('/home/lufla/workspace/brams-wind/Shapefiles/HIntegrada1', 'rios',color='w')
	#Fill continents
	if fillcont=='1':	
		m.fillcontinents(color='#BDA973')
	return

def doBasicMap(totx,toty,centlat,centlon):
	#Creating a map with data above
	m = Basemap(width=totx,height=toty,lon_0=centlon,lat_0=centlat, projection='lcc',resolution='h')
	m.drawcoastlines()
	lbs=[True, False, False, True]
	m.drawmeridians(range(-180,180,1),labels=lbs);
	m.drawparallels(range(-90,90,1),labels=lbs);
	return m

def readGrads(gradsf,m):
	ga = GaNum(Bin='grads',Window=False)
	try:
		fh = ga.open(gradsf+".ctl")
	except:
		print 'ERROR: grads file not found. Check if file '+gradsf+' (ctl and bin) is available'
		sys.exit(0)
	ts = ga.exp("tempk")
	nx=len(ts.grid.lon)
	ny=len(ts.grid.lat)
	x = linspace(0, m.urcrnrx, nx)
	y = linspace(0, m.urcrnry, ny)
	xa,ya = meshgrid(x,y)
	return ga,xa,ya

def plotTempkTopo(ga,xa,ya,m):
	ts = ga.exp("tempk")
	tempk=ts-273.0
	cs = m.contourf(xa,ya,tempk,30,linewidths=1.5,cmap=plt.cm.spring)
	cbar = m.colorbar(cs,location='bottom',pad="5%")
	unid=u"\U000000B0"+"C"
	cbar.set_label(unid)
	top=ga.exp("topo")
	topo=top
	ts = m.contour(xa,ya,topo,15,linewidths=0.5,colors="gray")
	plt.clabel(ts, inline=1, fontsize=5,fmt='%d',linestyles='dashed')
	plt.title("Temperature & Topo")
	return

def plotWinds(ga,xa,ya,m,plotarrow):
	uwind=ga.exp("u")
	vwind=ga.exp("v")
	wwind=ga.exp("w")
	u=uwind
	v=vwind
	w=wwind
	speed = np.sqrt(u*u + v*v)
	ww=m.contourf(xa,ya,w,30,linewidths=1.5,cmap=plt.cm.Wistia)
	wbar = m.colorbar(ww,location='right',pad="20%")
	unid="m/s"
	wbar.set_label("W (Shaded) "+unid)
	
	
	if plotarrow=='1':
		ws = m.quiver(xa,ya,u,v,speed,cmap=plt.cm.cool)
		wbar2 = m.colorbar(ws,location='bottom',pad="5%")
		wbar2.set_label("UV (vectors) "+unid)
	else:
		lw = 5*speed / speed.max()
		ws=m.streamplot(xa, ya, u, v,linewidth='1.0',color=speed,cmap=plt.cm.cool)
		wbar2 = m.colorbar(ws.lines,location='bottom',pad="5%")
		wbar2.set_label("UV (streamlines) "+unid)
	plt.title("Winds uv & w")
	return

def plotTurb(plotturbines,xt,yt,m,turb,nturb):
	if plotturbines=='1':
		i=0
		while i < nturb:
			x1,y1=m([xt[i]],[yt[i]])
			m.plot(x1,y1, marker="^", label='#'+turb[i])
			i=i+1
		plt.legend()

def readHdf5(param):
	filename = param+'.h5'
	try:
		file = h5py.File(filename, 'r')
	except:
		print 'ERROR: File not found. Check if file '+filename+' is available'
		sys.exit(0)
	
	#looking for name of windfarm and turbine (attributes)
	item = file["File info"]
	windfarm=item.attrs['Windfarm']
	turbine=item.attrs['#Turbine']
	
	# List all groups
	#print("Keys: %s" % file.keys())
	
	# Get the lat lon of turbine
	a_group_key = file.keys()[1]
	pgeo = list(file[a_group_key])
	#print 'Model lat: ',pgeo[0],', real lat: ',pgeo[2]
	#print 'Model lon: ',pgeo[1],', real lon: ',pgeo[3]
	
	#####################################
	#Getting hdf5Information from hdf5 file
	#####################################
	# Get the date/filehdf5Info
	a_group_key = file.keys()[0]
	hdf5Info = list(file[a_group_key])
	
	
	# Get the power
	a_group_key = file.keys()[2]
	power = list(file[a_group_key])
	
	# Transform datetime in datetime format of Python
	hhour=int(str(hdf5Info[3])[0:1])
	hmin=int(str(hdf5Info[3])[2:3])
	hsec=int(str(hdf5Info[3])[4:5])
	tini=datetime.datetime(hdf5Info[0],hdf5Info[1],hdf5Info[2],hhour,hmin,hsec)
	hhour=int(str(hdf5Info[9])[0:1])
	hmin=int(str(hdf5Info[9])[2:3])
	hsec=int(str(hdf5Info[9])[4:5])
	tend=datetime.datetime(hdf5Info[6],hdf5Info[7],hdf5Info[8],hhour,hmin,hsec)
	
	# Get the direction of wind
	a_group_key = file.keys()[3]
	dir = list(file[a_group_key])
	
	# Get the u wind 
	a_group_key = file.keys()[4]
	u = list(file[a_group_key])
	
	# Get the v wind 
	a_group_key = file.keys()[5]
	v = list(file[a_group_key])
	
	# Get the total wind 
	a_group_key = file.keys()[6]
	wind = list(file[a_group_key])
	
	t=np.arange(0,len(power),1)
	
	return power,tini,tend,u,v,wind,dir,t,windfarm,turbine,pgeo,hdf5Info

def writeCSV(param,power,tini,u,v,wind,dir,hdf5Info):
	fcvs=open(param+'.csv', 'wt')
	writer = csv.writer(fcvs,delimiter =";",)
	writer.writerow( ('Time', 'Power [W]','U wind [m/s]','v_wind [m/s]','Total Wind [m/s]','Wind dir [dg from N]') )
	for x in range(0,len(power)):
		datM=tini+datetime.timedelta(0,hdf5Info[5]*x)
		writer.writerow((datM.strftime("%Y/%m/%d %H:%M:%S"),power[x],u[x],v[x],wind[x],dir[x]))
	fcvs.close()
	return

def plotPower(t,power,hdf5Info,windfarm,turbine,pgeo,tini,tend):
	plt.plot(t,power,label='power')
	plt.xlabel('Time from start x '+'{:02}'.format(hdf5Info[5])+' Sec')
	plt.ylabel('Power [W]')
	plt.grid(True)
	plt.title(windfarm[0]+' '+'#'+turbine[0]+'{:12.6f}'.format(pgeo[0])+' '+'{:12.6f}'.format(pgeo[1])+'\n'+
		'Start: '+str(tini)+' End: '+str(tend))
	return

def plotuv(t,u,v,wind,hdf5Info):
	plt.plot(t,u,label='u')
	plt.plot(t,v,label='v')
	plt.plot(t,wind,label='total')
	plt.xlabel('Time from start x '+'{:02}'.format(hdf5Info[5])+' Sec')
	plt.ylabel('Velocity [m/s]')
	plt.grid(True)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=9, borderaxespad=0.)
	return

def plotWindRose(dir,wind,windfarm,turbine,pgeo,tini,tend):
	ax = WindroseAxes.from_ax()
	ax.bar(dir,wind,normed=True, opening=0.8, edgecolor='white')
	ax.set_legend()
	plt.title(windfarm[0]+' '+'#'+turbine[0]+'{:12.6f}'.format(pgeo[0])+' '+'{:12.6f}'.format(pgeo[1])+'\n'+
	          'Start: '+str(tini)+' End: '+str(tend))
	return

if __name__ == "__main__":
   main(sys.argv[1:])