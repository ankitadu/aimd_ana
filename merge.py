#!/usr/bin/env python 
#

import os
writedir = os.getcwd()
from decimal import *
#getcontext().prec = 5

#def skip_header():
#	if line_i1.startswith("#") or line_i1.startswith("!"):
#		print ("skipping header")
#		print (line_i1)
#		return True


i1 = open('data_equil.txt','r')	#filter template
i2 = open('temp3.txt','r')								#material to be filtered
content2 = i2.readlines()

o = open('temp.txt','w')
o.write('# headlines and comments\n')

attemptcounter = 0
found_hop_flag = 0

for line_i1 in i1:
	if found_hop_flag == 0 and attemptcounter > 0:
		print ("couldn't find matching frame or ihyd ")
	nstep, ihyd, dehop = line_i1.split()
	attemptcounter = 0
	found_hop_flag = 0
#	found_dist_flag = 0
	for line_i2 in content2:
		#if line_i2.startswith(frame+"\t"):#CAVE: space or tab?
			#attemptcounter += 1
#			print ("found frame "+frame+" hop attempt " +str(attemptcounter))
			#frame2,dEhop2,jres2,rhop2,oohop2,covb2,dh32,dphi2,theto2,nofs2,nnfs2,dein2,deac2,defs2,deenv2,dhb2,doo2,dthet2,devhb2,devoo2,devth2 = line_i2.split()
			#frame2,AIMD_ires,CHARMM_jres,dEhop2,rhop2 = line_i2.split()
			nframe, nhyd, angle, distanceOO, distanceHdum = line_i2.split()

			if (nframe == nstep):
				if (nhyd == ihyd):
					o.write(nstep+"\t"+ihyd+"\t"+dehop+"\t"+angle+"\t"+distanceOO+"\t"+distanceHdum+"\n")
					#found_hop_flag = 1	
				#elif ( round(Decimal(rhop),4) == round(Decimal(rhop2),4) ):
					#o.write(line_i2)
					#found_hop_flag = 1
				#elif ( round(Decimal(rhop),3) == round(Decimal(rhop2),3) ):
					#o.write(line_i2)
					#found_hop_flag = 1
				#elif ( round(Decimal(rhop),2) == round(Decimal(rhop2),2) ):
					#o.write(line_i2)
					#found_hop_flag = 1
			#else:
				#if attemptcounter > 3 :
					#print ("more than three hops to unrequested acceptor in frame "+frame)
					#exit()		

print ("filtering completed")
i2.close()
i1.close()
o.close()

