include(FortranCInterface)

FortranCInterface_VERIFY(CXX)

set(SOURCE 
    ${SOURCE}
    ${CMAKE_CURRENT_SOURCE_DIR}/itinmcm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/itinmcm.f90
    PARENT_SCOPE
)

set(HEADER
    ${HEADER}
    ${CMAKE_CURRENT_SOURCE_DIR}/itinmcm.hpp
)