find_package(Qt5Widgets)
find_package(Qt5PrintSupport)

if(Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND)
    message(STATUS "Found Qt5::Widgets and Qt5::PrintSupport")
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
    add_definitions(-DUSE_GUI ${Qt5Widgets_DEFINITIONS})
    qt5_wrap_ui(AUTORALLYRT_UI_HDRS autorallyrtwindow.ui)
    add_executable(AutoRallyRealTime autorallyrt.cpp
            ${CMAKE_SOURCE_DIR}/examples/qcustomplot/qcustomplot.cpp
            autorally_bfm.h
            autorally_costs.h
            autorallyrtwindow.h
            autorallyrtwindow.cpp
            ${AUTORALLYRT_UI_HDRS})
    target_include_directories(AutoRallyRealTime PUBLIC
            ${CMAKE_SOURCE_DIR}/include
            ${Qt5Widgets_INCLUDE_DIRS}
            ${Qt5PrintSupport_INCLUDE_DIRS}
            ${CMAKE_SOURCE_DIR}/examples)
    target_link_libraries(AutoRallyRealTime Qt5::Widgets Qt5::PrintSupport)
    target_compile_options(AutoRallyRealTime PRIVATE $<$<CONFIG:Release>:-ffast-math>)
    add_custom_command(TARGET AutoRallyRealTime POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy
                       ${CMAKE_CURRENT_SOURCE_DIR}/THETA.txt $<TARGET_FILE_DIR:AutoRallyRealTime>)
else()
    message(WARNING "Qt5 not found -- AutoRally simulation not built")
endif()
