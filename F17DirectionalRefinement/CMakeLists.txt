
add_pz_target(
  NAME F17DirRef
  SOURCES f17.cpp
  FILES
  yf17.xyz
  yf17.tri
  yf17.tet)
#in this directory there are also the original mesh files
#zipped in a file. the script truncateMesh.py could be used
#for creating a smaller bounding box