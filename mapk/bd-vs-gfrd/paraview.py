base_directory = '/home/miedema/code/amolf/'
egfrd_directory = base_directory + 'egfrd/'

#simulation_data_directory = base_directory + 'simulations/mapk/bd-vs-gfrd/data/bd-1800_cpu_s-1.0_kk_vs_p/'
#simulation_data_directory = base_directory + 'simulations/mapk/bd-vs-gfrd/data/gfrd-1800_cpu_s-1.0_kk_vs_p/'


import sys
sys.path.append(egfrd_directory)

from visualization import pipeline


# Create pipeline.
p = pipeline.Pipeline(egfrd_directory, simulation_data_directory)
p.particle_radius_scale_factor = 10
p.resolution = 8
p.rebuild()

