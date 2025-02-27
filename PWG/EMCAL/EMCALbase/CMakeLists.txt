# **************************************************************************
# * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: The ALICE Off-line Project.                                    *
# * Contributors are mentioned in the code where appropriate.              *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************/

#Module
set(MODULE PWGEMCALbase)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWG/EMCAL/EMCALbase)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/TENDER/Tender
                    ${AliPhysics_SOURCE_DIR}/TENDER/TenderSupplies
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools/yaml-cpp/include
  )

# Sources - alphabetical order
set(SRCS
  AliEMCALRecoUtils.cxx
  AliAnalysisTaskEmcal.cxx
  AliAnalysisTaskEmcalLight.cxx
  AliClusterContainer.cxx
  AliEmcalContainer.cxx
  AliEmcalContainerUtils.cxx
  AliEmcalDownscaleFactorsOCDB.cxx
  AliEmcalCutBase.cxx
  AliEmcalVCutsWrapper.cxx
  AliEmcalAODFilterBitCuts.cxx
  AliEmcalAODHybridTrackCuts.cxx
  AliEmcalAODTPCOnlyTrackCuts.cxx
  AliEmcalContainerUtils.cxx
  AliEmcalESDTrackCutsGenerator.cxx
  AliEmcalESDHybridTrackCuts.cxx
  AliEmcalESDtrackCutsWrapper.cxx
  AliEmcalParticle.cxx
  AliEmcalPhysicsSelection.cxx
  AliEmcalPythiaInfo.cxx
  AliEmcalTrackSelResultPtr.cxx
  AliEmcalTrackSelResultCombined.cxx
  AliEmcalTrackSelResultHybrid.cxx
  AliEmcalTrackSelection.cxx
  AliEmcalTrackSelectionESD.cxx
  AliEmcalTrackSelectionAOD.cxx
  AliParticleContainer.cxx
  AliPicoTrack.cxx
  AliMCParticleContainer.cxx
  AliTrackContainer.cxx
  AliEmcalList.cxx
  AliAnalysisTaskEmcalEmbeddingHelper.cxx
  AliEmcalEmbeddingQA.cxx
  )


if(RooUnfold_FOUND)
    get_target_property(ROOUNFOLD_INCLUDE_DIR RooUnfold::RooUnfold INTERFACE_INCLUDE_DIRECTORIES)
    include_directories(${ROOUNFOLD_INCLUDE_DIR})
    add_definitions("-DWITH_ROOUNFOLD")
endif(RooUnfold_FOUND)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
set(HDRS
  "${HDRS}"
  AliEmcalIterableContainer.h
  AliEmcalContainerIndexMap.h
  AliEmcalStringView.h
  )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} ${HDRS} ${MODULE}LinkDef.h G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSIS ANALYSISalice AOD OADB CDB EMCALrec EMCALUtils ESD PWGTools STEER STEERBase Tender)
if(RooUnfold_FOUND)
    get_target_property(ROOUNFOLD_LIBRARY RooUnfold::RooUnfold IMPORTED_LOCATION)
    set(LIBDEPS ${LIBDEPS} ${ROOUNFOLD_LIBRARY})
endif(RooUnfold_FOUND)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

# Unit tests

add_test(func_PWGEMCALbase_AliEmcalTrackSelResultPtr
    env
    LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
    DYLD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{DYLD_LIBRARY_PATH}
    ROOT_HIST=0
    root -n -l -b -q "${CMAKE_INSTALL_PREFIX}/PWG/EMCAL/macros/TestAliEmcalTrackSelResultPtr.C)")

add_test(func_PWGEMCALbase_AliEmcalAODHybridTrackCuts
    env
    LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
    DYLD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{DYLD_LIBRARY_PATH}
    ROOT_HIST=0
    root -n -l -b -q "${CMAKE_INSTALL_PREFIX}/PWG/EMCAL/macros/TestAliEmcalAODHybridTrackCuts.C)")

add_test(func_PWGEMCALbase_AliEmcalTrackSelectionAOD
    env
    LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
    DYLD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:$ENV{DYLD_LIBRARY_PATH}
    ROOT_HIST=0
    root -n -l -b -q "${CMAKE_INSTALL_PREFIX}/PWG/EMCAL/macros/TestAliEmcalTrackSelectionAOD.C)")
    
