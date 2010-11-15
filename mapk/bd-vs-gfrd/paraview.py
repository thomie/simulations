base_directory = '/home/miedema/code/amolf/'
paraview_scripts_directory = base_directory + 'egfrd/visualization/'
#simulation_data_directory = base_directory + 'simulations/bd-vs-gfrd/data/bd-1800_cpu_s-1.0_kk_vs_p/'
simulation_data_directory = base_directory + 'simulations/bd-vs-gfrd/data/gfrd-1800_cpu_s-1.0_kk_vs_p/'

execfile(paraview_scripts_directory + 'pipeline.py')
