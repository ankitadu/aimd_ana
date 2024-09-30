# aimd_ana
These codes analyze specific ab initial trajectories and compare them to outputs from Classical MD trajectories 
'WaterHydroxide27Images_wrapped_Center.f90' : Will read through Qchem trajectory data and then write cartesian coordinates of the molecules of interest for visualization and structural analysis. 

'3_ana_f-dip_backgeo_wrapper.py' : is able to create CHARMM MD inputs flies using the previous coordinates and then run a short simulation. This calculates the energetic and structural difference between Qchem vs classical MD- The ensemble average of which is our target function. 

'hop_count.py' : is able to calculate the number of proton hops by counting the number of times a molecular ID changes. Proton hops occur in burst so we don't want to count continuous back and forth hop (rattling) this code also stores previous IDs to eliminate such rattling from the actual proton hop count. 

'merge.py' : simple python code to clean up the data and make sure there aren't any repeating data points

'derivativescopy_norm.f90' : is able to take in all the structure/confirmation and energetic data produced previously and calculates multi-parameter gradients to minimize our target function using steepest descent. We use a simple bash executable to run MD on CHARMM iteratively until a minima is reached. 
