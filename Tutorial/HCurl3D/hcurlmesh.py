import os
import gmsh
import sys

gmsh.initialize()

gmsh.model.add("hcurlmesh")

gmsh.logger.start()

# volume identifier
vol_tag = 1
# boundary identifier
bnd_tag = 2
# element size
el_size = 0.3
# x, y, z, dx, dy, dz
gmsh.model.occ.addBox(-1, -1, -1, 2, 2, 2, vol_tag)
gmsh.model.occ.synchronize()
bndlist = [t for _, t in gmsh.model.get_boundary([(3,vol_tag)], oriented=False)]
gmsh.model.add_physical_group(3,[vol_tag], vol_tag, "vol")
gmsh.model.add_physical_group(2,bndlist, bnd_tag, "bnd")

# Assign a mesh size to all the points:
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), el_size)

gmsh.model.mesh.generate(3)

if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    gmsh.write("hcurlmesh.msh")

# Additional examples created with the OpenCASCADE geometry kernel are available
# in `t18.py', `t19.py' and `t20.py', as well as in the `examples/api'
# directory.

# Inspect the log:
log = gmsh.logger.get()
print("Logger has recorded " + str(len(log)) + " lines")
gmsh.logger.stop()

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

