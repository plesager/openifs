# (C) Copyright 1996-2019 ECMWF.

get_filename_component(IFS_EXECUTABLE_NAME ${IFS_EXECUTABLE} NAME)

if( TARGET ${IFS_EXECUTABLE_NAME} )
  list( APPEND TEST_DEPENDS ${IFS_EXECUTABLE_NAME} )
endif()

if( IFS_EXECUTABLE_NAME STREQUAL ifsMASTER.DP )
  set(DOUBLE_PRECISION ON)
endif()

set(APATH_255l_2 /ec/res4/hpcperm/nld5163/openifs/oifs_data/48r1/climate.v020/255l_2) 

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/ifsdata ${CMAKE_CURRENT_BINARY_DIR}/ifsdata )

add_subdirectory( common )

#################################################################
# Helper function
function( add_ifs_test TARGET )

    ### (1) COMMON for all tests 
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/${TARGET} )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-run ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-run )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-debug ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-debug )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-clean ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-clean )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-check-bitid ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-check-bitid )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-check-tracers ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-check-tracers )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/ifs-grep-norms.pl ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs-grep-norms.pl )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${PROJECT_SOURCE_DIR}/bin/set_launcher.bash ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/set_launcher.bash )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${SHARE_BINARY_DIR}/ifs_env ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifs_env )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/inidata/ICMGGham7INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGham7INIT )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/inidata/ICMGGham7INIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGham7INIUA )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/inidata/ICMSHham7INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHham7INIT )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/inidata/ICMCLham7INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLham7INIT )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/inidata/ICMCLham7INIT_COMPO ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLham7INIT_COMPO )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/m7/lut_optical_properties_lw_M7.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/lut_optical_properties_lw_M7.nc)
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/m7/lut_optical_properties_M7.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/lut_optical_properties_M7.nc)

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/m7/lookup_table.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/lookup_table.nc)
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/m7/refractive_indices_hdfstyle.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/refractive_indices_hdfstyle.nc)
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/m7/namelist.echam ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/namelist.echam)
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/m7/parnuc.15H2SO4.A0.total.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/parnuc.15H2SO4.nc)
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/m7/lut_kappa.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/lut_kappa.nc)

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/OMI.data.extraterrest ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/OMI.data.extraterrest )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/aerosol_reduce.dat ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/aerosol_reduce.dat )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/tropo_look_up_cbmhybrid.dat ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/tropo_look_up_cbmhybrid.dat )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/uars_ratio.txt ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/uars_ratio.txt )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/haloe_ch4clim.dat ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/haloe_ch4clim.dat )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${APATH_255l_2} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/255l_2 )

    # execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    #  ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data/255l_2 ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/255l_2 )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/ifsdata ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifsdata )

    if( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/params" )
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/params ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/params )
    endif()

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/setup ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/setup )

    if( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/postprocessing" )
      execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}/postprocessing ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/postprocessing )
    endif()

    ecbuild_add_test( TARGET ifs_t255_${TARGET}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}
              COMMAND "./ifs-run"
              DEPENDS ${TEST_DEPENDS}
            )
    set_property(TEST ifs_t255_${TARGET} APPEND PROPERTY LABELS t255)

endfunction()

#################################################################
# add tests

add_ifs_test( test_compo_fc_m7 )

# Test which can only be run if OOPS-IFS has been compiled:
# if( ifs_HAVE_OOPS )
#  add_ifs_test( test_compo_fc_dualconf )
# endif()


if( ENABLE_TESTS )

  set_property(TEST ifs_t255_test_compo_fc_m7 APPEND PROPERTY LABELS forecast mpi threading composition m7)

endif()
