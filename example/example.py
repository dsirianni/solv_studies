import sys
sys.path.insert(1,'/Users/bbakr/Desktop/gits/solv_studies/')

from geom import geom


#Example
geometry = geom('geom.xyz', solvent = 'HOH')
geometry.write_input(file_name = 'hello_world.in', solv_mon = 'B')

##to put solvent in C by default where all solvent is water for list of geometries 
#for geometry in geometry_list:
#    this_geom = geom(geometry, solvent = 'HOH')
#    this_geom.write_ipnut(file_name = geometry.split('.xyz')[0]+'.in', solv_mon = 'C')
