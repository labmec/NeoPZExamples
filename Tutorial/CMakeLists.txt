set(tutorial_sources "")

add_subdirectory(Poisson2D)
add_subdirectory(Projections)


#enable folder structure for the tutorials (needed in some IDEs)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${tutorial_sources})