
if(USE_THERMOFUN MATCHES ON)
  find_package(ThermoFun REQUIRED)
  if(NOT ThermoFun_FOUND)
    message(FATAL_ERROR "ThermoFun library not found")
  else()
    message(STATUS "Found ThermoFun v${ThermoFun_VERSION}")
  endif()

  find_package(ChemicalFun REQUIRED)
  if(NOT ChemicalFun_FOUND)
    message(FATAL_ERROR "ChemicalFun library not found")
  else()
    message(STATUS "Found ChemicalFun v${ChemicalFun_VERSION}")
  endif()

endif()

if(USE_SPDLOG_PRECOMPILED)
 if(NOT TARGET spdlog::spdlog)
    find_package(spdlog CONFIG REQUIRED)
    if(NOT spdlog_FOUND)
       message(FATAL_ERROR "spdlog library not found")
    else()
       message(STATUS "Found spdlog v${spdlog_VERSION}")
    endif()
  endif()
endif()

find_package(GEMS3K REQUIRED)
if(NOT GEMS3K_FOUND)
   message(FATAL_ERROR "GEMS3K library not found")
else()
   message(STATUS "Found GEMS3K v${GEMS3K_VERSION}")
endif()
