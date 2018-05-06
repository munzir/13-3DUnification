find_package(Qt5Widgets)
find_package(Qt5PrintSupport)

if(Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND)
    message(STATUS "Found Qt5::Widgets and Qt5::PrintSupport")
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
    add_definitions(-DUSE_GUI ${Qt5Widgets_DEFINITIONS})
    qt5_wrap_ui(KRANG3D_UI_HDRS krang3dwindow.ui)
    add_executable(Krang3D krang3d.cpp
            ${CMAKE_SOURCE_DIR}/examples/qcustomplot/qcustomplot.cpp
            krang3dwindow.h
            krang3dwindow.cpp
            ${KRANG3D_UI_HDRS})
    target_include_directories(Krang3D PUBLIC ${Qt5Widgets_INCLUDE_DIRS}
            ${Qt5PrintSupport_INCLUDE_DIRS}
            ${CMAKE_SOURCE_DIR}/examples)
    target_link_libraries(Krang3D Qt5::Widgets Qt5::PrintSupport)
else()
    message(WARNING "Qt5 not found -- Krang3D simulation not built")
endif()
