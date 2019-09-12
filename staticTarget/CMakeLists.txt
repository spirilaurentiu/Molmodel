## This whole directory exists just so I could define these extra 
## preprocessor values.
ADD_DEFINITIONS(
  -DSimTK_MOLMODEL_BUILDING_STATIC_LIBRARY 
  -DSimTK_USE_STATIC_LIBRARIES)

IF(BUILD_UNVERSIONED_LIBRARIES)

    # This is the list of dependencies so we know when to rebuild.
    ADD_LIBRARY(${STATIC_TARGET} STATIC 
                ${SOURCE_FILES} 
                ${SOURCE_INCLUDE_FILES} 
                ${API_ABS_INCLUDE_FILES})
    
    TARGET_LINK_LIBRARIES(${STATIC_TARGET}
	${Simbody_STATIC_LIBRARIES}
        ${ADDITIONAL_LINK_LIBRARIES})
    
    SET_TARGET_PROPERTIES(${STATIC_TARGET} PROPERTIES
    	PROJECT_LABEL "Library - ${STATIC_TARGET}")
    
    # install library; on Windows .dll goes in the bin directory.
    INSTALL_TARGETS(/lib${LIB64} 
                    RUNTIME_DIRECTORY /bin
                    ${STATIC_TARGET})

ENDIF(BUILD_UNVERSIONED_LIBRARIES)



IF(BUILD_VERSIONED_LIBRARIES)

    # This is the list of dependencies so we know when to rebuild.
    ADD_LIBRARY(${STATIC_TARGET_VN} STATIC 
                ${SOURCE_FILES} 
                ${SOURCE_INCLUDE_FILES} 
                ${API_ABS_INCLUDE_FILES})
    
    TARGET_LINK_LIBRARIES(${STATIC_TARGET_VN}
	${Simbody_STATIC_LIBRARIES_VN}
        ${ADDITIONAL_LINK_LIBRARIES})
    
    SET_TARGET_PROPERTIES(${STATIC_TARGET_VN} PROPERTIES
    	PROJECT_LABEL "Library - ${STATIC_TARGET_VN}")
    
    # install library; on Windows .dll goes in the bin directory.
    INSTALL_TARGETS(/lib${LIB64} 
                    RUNTIME_DIRECTORY /bin
                    ${STATIC_TARGET_VN})

ENDIF(BUILD_VERSIONED_LIBRARIES)
