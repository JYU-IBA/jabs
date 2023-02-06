set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(../)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/..)


add_executable(des_test
        des_test.c
        ../sample.c ../brick.c ../ion.c ../simulation.c ../reaction.c
        ../spectrum.c ../fit.c ../fit_params.c ../rotate.c ../detector.c ../jabs.c
        ../roughness.c  ../script.c  ../generic.c ../message.c ../aperture.c
        ../geostragg.c  ../script_command.c ../script_session.c ../script_file.c
        ../calibration.c ../prob_dist.c ../nuclear_stopping.c ../stop.c ../des.c
        ../simulation_workspace.c ../sim_reaction.c ../sim_calc_params.c
        "$<$<BOOL:${JABS_PLUGINS}>:../plugin.c>"
        ../idfparse.c ../idf2jbs.c ../idfelementparsers.c ../options.c
)

target_link_libraries(des_test
    PRIVATE jibal
    PRIVATE GSL::gsl
    PRIVATE LibXml2::LibXml2
    PRIVATE "$<$<BOOL:${UNIX}>:m>"
    gitwatcher
)
