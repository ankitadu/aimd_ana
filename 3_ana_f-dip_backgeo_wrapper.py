#!/usr/bin/env python 
#

import os
import time

time_to_wait = 5 # seconds
time_counter = 0

#choose one replacement term (lists start at slice number[0] )
models = ["sv","pol"]
model = models[0]

#resids = [17, 27, 29]#for testing
resids = range(1,32) #make list from 1 to 31, 32 is donor!

outdir = "out_dipole_rbackgeo/"
for resid in resids:
	with open(outdir+"hydr"+str(resid)+".ts", "w") as o:#write header
		o.write("#tilt angle betw mol dipole + H3O's O, pre- and post-hop, distance O...closest donor H (rhop)\n")
		o.write("#frame\tO...H dip_magn thet_oh2 f_thet_oh2 dihH_thet f_rhop H3H3 thet_H3H3O\n")

with open ("template_addons_f-dipa_rbackgeo_H.inp","r") as i0:
	template = i0.readlines()

for resid in resids:
	print ("analyzing residue nr "+str(resid))
	#framecounter sets a limit for test purposes
	framecount = 0
	while framecount < 15299 :#this script takes ~13 hours to run on all 15300 frames 15299
		framecount += 1
		o = open(outdir+"hydr"+str(resid)+".ts", "a") #write time series, append mode
		o.write(str(framecount)+"\t")
		o.close() #need to close so that result can be appended properly to closed file later

		#write new input file for each frame
		infile = model+"_ana_f-dip_rbackgeo_H.inp"
		cast = open(infile,"w")

		terms_to_replace = ["MODEL","xRESID","FRAMECOUNT"]
		replacements = [model,str(resid),str(framecount)]

		for t_line in template:
			cast_line = t_line
			for term in terms_to_replace:
				idx = terms_to_replace.index(term)
				cast_line = cast_line.replace(term,replacements[idx])
			cast.write(cast_line)
		cast.close()
		#process each individual frame with CHARMM
		outfile = model+'_frame_'+str(framecount)+'.out'
		os.system("charmm < "+model+"_ana_f-dip_rbackgeo_H.inp > "+outdir+outfile)

		anafile_path = outdir+model+"_ana_f-dipole_too_rbackgeo_h_hydr"+str(resid)+"_frame_"+str(framecount)+".dat"
		while not os.path.exists(anafile_path): #careful with capital letters, CHARMM can't write them in file names apparently
   			time.sleep(1)
			print("waiting")
   			time_counter += 1
   			if time_counter > time_to_wait:
				print ("analysis of frame "+str(framecount)+" takes too long, exiting")
				exit()

		os.system ("tail -1 "+anafile_path+" >> "+outdir+"hydr"+str(resid)+".ts") #append output to collection immediately
		os.system ("rm "+anafile_path)		 #clean up immediately, so I don't get a gazillion files
		os.system("rm "+outdir+outfile)

print( "done with "+str(len(resids)*framecount) +" readings")
