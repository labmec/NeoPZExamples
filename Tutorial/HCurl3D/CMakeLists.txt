add_pz_target(
    NAME HCurl3D
    SOURCES HCurl3D.cpp TPZMatHCurl3D.cpp
    FILES hcurlmesh.msh)

list(APPEND tutorial_sources HCurl3D.cpp TPZMatHCurl3D.cpp TPZMatHCurl3D.h)