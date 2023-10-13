if (COMPILE_ITSLAV)
    add_subdirectory(itslav)
endif()

if (COMPILE_ITINMCM)
    add_subdirectory(itinmcm)
endif()

set(SOURCE ${SOURCE} PARENT_SCOPE)
set(HEADER ${HEADER} PARENT_SCOPE)