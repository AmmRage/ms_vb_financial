# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

# 项目信息
project (ms_vb_financial)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 添加 math 子目录
add_subdirectory(gtest)

# 指定生成目标
add_executable(ms_vb_financial test_main.cc)



include_directories(gtest/googletest/gtest)
include_directories(src)
# 添加链接库
#target_link_libraries(ms_vb_financial 
#	${CMAKE_CURRENT_SOURCE_DIR}/lib/gtestd.lib
#)
target_link_libraries(ms_vb_financial gtest)

if (WIN32)
	set_property(TARGET ms_vb_financial PROPERTY
	  MSVC_RUNTIME_LIBRARY "MultiThreadedDebug$<$<CONFIG:Debug>:Debug>")
endif (WIN32)