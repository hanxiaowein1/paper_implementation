cmake_minimum_required(VERSION 3.20.0)

set(CMAKE_CXX_STANDARD 20)
set(LIBIGL_WITH_CGAL ON)
set(CMAKE_TOOLCHAIN_FILE "D:\\Library\\vcpkg\\scripts\\buildsystems\\vcpkg.cmake")
# add_compile_options must added before add_executable and add_library...
add_compile_options(/bigobj)

project(low_poly_meshing)
find_package(GTest CONFIG REQUIRED)
find_package(Boost COMPONENTS functional REQUIRED)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} Boost::functional)
target_link_libraries(${PROJECT_NAME} GTest::gtest GTest::gmock)

# for unit test
enable_testing()
find_package(CGAL REQUIRED)
set(IGL_STATIC_LIBRARY_DIR "D:\\Library\\libigl\\build\\lib\\Debug")
add_executable(tests unit_test.cpp feature_point_mc33.cpp feature_point_mc33_tables.cpp)
target_link_libraries(tests GTest::gtest GTest::gmock)
target_link_libraries(tests CGAL::CGAL)
target_link_libraries(tests Boost::functional)
link_directories(IGL_STATIC_LIBRARY_DIR)
add_test(AllTestsInMain tests)

include_directories("D:\\Library\\libigl\\build\\_deps\\eigen-src")
include_directories("D:\\Library\\libigl\\include")
include_directories("D:\\Library\\libigl\\build\\_deps\\glfw-src\\include")
include_directories("D:\\Library\\libigl\\build\\_deps\\glad-src\\include")
