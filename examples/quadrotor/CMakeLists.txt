find_package(Qt5Widgets)
find_package(Qt5PrintSupport)

if(Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND)
    message(STATUS "Found Qt5::Widgets and Qt5::PrintSupport")
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
    add_definitions(-DUSE_GUI ${Qt5Widgets_DEFINITIONS})
    qt5_wrap_ui(QUADROTOR_UI_HDRS quadrotorwindow.ui)
    add_executable(Quadrotor quadrotor.cpp quadrotor.h quadrotor_costs.h quadrotor_plant.h
            ${CMAKE_SOURCE_DIR}/examples/qcustomplot/qcustomplot.cpp
            quadrotorwindow.h
            quadrotorwindow.cpp
            ${QUADROTOR_UI_HDRS})
    target_include_directories(Quadrotor PUBLIC ${Qt5Widgets_INCLUDE_DIRS}
            ${Qt5PrintSupport_INCLUDE_DIRS}
            ${CMAKE_SOURCE_DIR}/examples)
    target_link_libraries(Quadrotor Qt5::Widgets Qt5::PrintSupport)
else()
    message(WARNING "Qt5 not found -- Quadrotor simulation not built")
endif()
