add_pz_target(
    NAME HCurlProjection
    SOURCES HCurlProjection.cpp)
add_pz_target(
    NAME HDivKernelProjection
    SOURCES HDivKernelProjection.cpp)

list(APPEND tutorial_sources HCurlProjection.cpp HDivKernelProjection.cpp)