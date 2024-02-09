from numpy import *
# for date functions
import datetime
# For system comand
import sys, os, glob, re

INIDATE=int(os.environ['INIDATE'])
###DIR_DATA="/home/oxana/pub2/ARCRES/FORECAST/MOLOCH/METEOGRAMS"
DIR_DATA=os.environ['DIR_METEOGRAMS']

###DIR_CURRENT=os.getcwd()
DIR_OUTPUT=DIR_DATA

AINIDATE=("%10.10i")%(INIDATE)

#os.chdir(DIR_DATA)

FILE_LIST=glob.glob(DIR_DATA+"/meteogram_"+AINIDATE+"_*.dat")
FILE_LIST.sort()

# Definition of NDATA:
FILEIN=open(FILE_LIST[0],'r')
LINE=FILEIN.readline()
LINE=FILEIN.readline()
FILEIN.close()
WORDS=LINE.split()
NDATA=len(WORDS)-4 # 4 first words are parameters ot the point

NFILE=len(FILE_LIST); NINST=NFILE
DATE0=[]; HFORECAST=[]
POINT_NAME=[]; POINT_LAT=[]; POINT_LON=[]; POINT_ALT=[]  
DATA=empty ( (1000, NINST, NDATA), dtype=float)
DATA1=empty ( (NDATA), dtype=float)
FILENAME_OUTPUT=[]

for IFILE in range(0,NFILE) :
  FILEIN=open(FILE_LIST[IFILE],'r')
  LINES=FILEIN.readlines()
  NLINES=len(LINES)
  LINE=LINES[0]
  WORDS=LINE.split()
  ####splitter = re.compile(' +'); WORDS=splitter.split(LINE)
  DATE0.append(datetime.datetime(int(WORDS[0]),int(WORDS[1]),int(WORDS[2]),int(WORDS[3]),0))
  HFORECAST.append(int(WORDS[4])*24+int(WORDS[5]))

  for ILINE in range(1,NLINES) :
    LINE=LINES[ILINE]
    WORDS=LINE.split()
    NAME=WORDS[0]
    LAT=float(WORDS[1])
    LON=float(WORDS[2])
    ALT=int(WORDS[3])
    for I in range (0,NDATA):
      DATA1[I]=float(WORDS[I+4])
    STRTIME=DATE0[IFILE].strftime("%Y%m%d%H")
    FILENAME_OUTPUT.append(("%s"*3+"%3.3i%s") % ("meteogram_", STRTIME, "_point_", ILINE, ".dat"))
   
    if IFILE == 0 : # first forecast instant
      
      POINT_NAME.append(NAME)
      POINT_LAT.append(LAT) 
      POINT_LON.append(LON) 
      POINT_ALT.append(ALT) 
      try : s=os.remove(FILENAME_OUTPUT[ILINE-1])
      except : pass
      FILEOUT=open(FILENAME_OUTPUT[ILINE-1],'w') 
      FIRST_LINE=("%s"+" %7.2f"*2+" %4i %s")%(NAME, LAT, LON, ALT, STRTIME)
      FILEOUT.write(FIRST_LINE+"\n")
      FILEOUT.close()
      NPOINT=len(POINT_NAME)

    else : # no first forecast instant

      FLAG=0
      for I in range(0,NPOINT):
        if NAME==POINT_NAME[I] and LAT==POINT_LAT[I] and LON==POINT_LON[I] and ALT==POINT_ALT[I] : FLAG=1; break
      if FLAG == 0 :
        print "MeteoStation not found in the list => pass"; break
      IPOINT=I 

    FILEOUT=open(FILENAME_OUTPUT[ILINE-1],'a') 
    FILEOUT.write(("%3.3i  "+"%7.1f"*NDATA+"\n")  % ( (HFORECAST[IFILE],) + tuple(DATA1[0:NDATA]) ))
    FILEOUT.close()

  FILEIN.close()
