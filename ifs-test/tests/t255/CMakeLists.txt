# (C) Copyright 1996-2019 ECMWF.

get_filename_component(IFS_EXECUTABLE_NAME ${IFS_EXECUTABLE} NAME)

if( TARGET ${IFS_EXECUTABLE_NAME} )
  list( APPEND TEST_DEPENDS ${IFS_EXECUTABLE_NAME} )
endif()

if( IFS_EXECUTABLE_NAME STREQUAL ifsMASTER.DP )
  set(DOUBLE_PRECISION ON)
endif()

# --- Paths for T255 datasets ---

# Determine inidata path (look in cache first, then source dir if extracted there)
find_path( T255_INIDATA_PATH
  NAMES ICMGGham7INIT
  PATHS
    ${DATA_CACHE_DIR}/compo_t255/ifs-test/tests/t255/inidata
    ${DATA_CACHE_DIR}/ifs-test/tests/t255/inidata
    ${DATA_CACHE_DIR}/t255/inidata
    ${CMAKE_CURRENT_SOURCE_DIR}/inidata
  NO_DEFAULT_PATH
)
if( NOT T255_INIDATA_PATH )
  set( T255_INIDATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/inidata" )
endif()

# Determine compo_data path
find_path( T255_COMPO_DATA_PATH
  NAMES aerosol_reduce.dat
  PATHS
    ${DATA_CACHE_DIR}/compo_t255/ifs-test/tests/compo_data
    ${DATA_CACHE_DIR}/ifs-test/tests/compo_data
    ${COMPO_DATA_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/../compo_data
  NO_DEFAULT_PATH
)
if( NOT T255_COMPO_DATA_PATH )
  set( T255_COMPO_DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../compo_data" )
endif()

# Determine static ifsdata path
find_path( T255_STATIC_DATA_PATH
  NAMES rtable_2021
  PATHS
    ${DATA_CACHE_DIR}/compo_t255/ifs-test/tests/t255/ifsdata
    ${DATA_CACHE_DIR}/ifs-test/tests/t255/ifsdata
    ${STATIC_DATA_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/ifsdata
  NO_DEFAULT_PATH
)
if( NOT T255_STATIC_DATA_PATH )
  set( T255_STATIC_DATA_PATH "${STATIC_DATA_PATH}" )
endif()

# Determine 255l_2 climate data directory
if( NOT APATH_255l_2 )
  if( DEFINED ENV{APATH_255l_2} )
    set( APATH_255l_2 "$ENV{APATH_255l_2}" )
  endif()
endif()

if( NOT APATH_255l_2 )
  set( _255l_2_search_paths
    "$ENV{OIFS_DATA_DIR}/$ENV{OIFS_CLIMATE}/255l_2"
    "$ENV{OIFS_DATA_DIR}/climate.v020/255l_2"
    "$ENV{OIFS_DATA_DIR}/48r1/climate.v020/255l_2"
    "$ENV{OIFS_DATA_DIR}/255l_2"
    "$ENV{OIFS_HOME}/openifs-data/climate.v020/255l_2"
    "$ENV{OIFS_HOME}/openifs-data/255l_2"
    "${PROJECT_SOURCE_DIR}/openifs-data/climate.v020/255l_2"
    "${PROJECT_SOURCE_DIR}/openifs-data/255l_2"
    "${DATA_CACHE_DIR}/climate.v020/255l_2"
    "${DATA_CACHE_DIR}/255l_2"
    "/hpcperm/nm6/oifs48r1-data/48r1/climate.v020/255l_2"
    "/ec/res4/hpcperm/nld5163/openifs/oifs_data/48r1/climate.v020/255l_2"
    "/projappl/project_465000527/models/ifs/ecpds/ifs-test-inidata/48r1/climate.v020/255l_2"
    "/gpfs/projects/bsc32/models/openifs/ifs-test-inidata/48r1/climate.v020/255l_2"
  )
  foreach( _path IN LISTS _255l_2_search_paths )
    if( EXISTS "${_path}" )
      set( APATH_255l_2 "${_path}" )
      break()
    endif()
  endforeach()
endif()

set( APATH_255l_2 "${APATH_255l_2}" CACHE PATH "Path to 255l_2 climate data directory" )

if( NOT EXISTS "${APATH_255l_2}" )
  ecbuild_warn("APATH_255l_2 was set to '${APATH_255l_2}', but the directory does not exist. Specify -DAPATH_255l_2=<path> or set OIFS_DATA_DIR in your environment.")
else()
  ecbuild_info("Using 255l_2 climate data from: ${APATH_255l_2}")
endif()

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${T255_STATIC_DATA_PATH} ${CMAKE_CURRENT_BINARY_DIR}/ifsdata )

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
      ${T255_INIDATA_PATH}/ICMGGham7INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGham7INIT )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_INIDATA_PATH}/ICMGGham7INIUA ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMGGham7INIUA )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_INIDATA_PATH}/ICMSHham7INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMSHham7INIT )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_INIDATA_PATH}/ICMCLham7INIT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLham7INIT )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_INIDATA_PATH}/ICMCLham7INIT_COMPO ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ICMCLham7INIT_COMPO )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/m7/lut_optical_properties_lw_M7.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/lut_optical_properties_lw_M7.nc)
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/m7/lut_optical_properties_M7.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/lut_optical_properties_M7.nc)

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/m7/lookup_table.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/lookup_table.nc)
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/m7/refractive_indices_hdfstyle.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/refractive_indices_hdfstyle.nc)
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/m7/namelist.echam ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/namelist.echam)
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/m7/parnuc.15H2SO4.A0.total.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/parnuc.15H2SO4.nc)
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/m7/lut_kappa.nc ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/lut_kappa.nc)

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/OMI.data.extraterrest ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/OMI.data.extraterrest )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/aerosol_reduce.dat ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/aerosol_reduce.dat )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/tropo_look_up_cbmhybrid.dat ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/tropo_look_up_cbmhybrid.dat )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/uars_ratio.txt ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/uars_ratio.txt )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/haloe_ch4clim.dat ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/haloe_ch4clim.dat )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/255l_2/L91/aeroclim.grib ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/aeroclim.grib )
    
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${APATH_255l_2} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/255l_2 )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_COMPO_DATA_PATH}/composition ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/composition )

    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${T255_STATIC_DATA_PATH} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}/ifsdata )

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
