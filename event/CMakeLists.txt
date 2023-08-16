
# Declare processing module
module( 
    NAME event 
    EXTERNAL_DEPENDENCIES ROOT LCIO 
)
target_include_directories(
  event PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>"
)
root_generate_dictionary(
  EventDict 
  EventDef.h 
  MODULE ${PROJECT_NAME} 
  LINKDEF ${event_INCLUDE_DIR}/EventLinkDef.h
)

# install ROOT pcm file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libevent_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libevent.rootmap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY ${event_INCLUDE_DIR} DESTINATION ${CMAKE_INSTALL_PREFIX})
