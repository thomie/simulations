
base_directory = '/home/miedema/code/amolf/'
egfrd_directory = base_directory + 'egfrd/'

#simulation_data_directory = base_directory + 'simulations/example/planar-surface_6-particles_200-steps/'
#simulation_data_directory = base_directory + 'simulations/example/12-particles_200-steps/'
#simulation_data_directory = base_directory + 'simulations/example/data/run'


simulation_data_directory = base_directory + 'simulations/example/self-evaluation-biochemical-networks-2005-2010/data/'


import sys
sys.path.append('/usr/lib/paraview')
sys.path.append(egfrd_directory)

from visualization import pipeline


# Create pipeline.
p = pipeline.Pipeline(egfrd_directory, simulation_data_directory)

p.particle_radius_scale_factor = 1
p.helix_radius_scale_factor = 1
p.cylinder_radius_scale_factor = 1
p.rebuild()
