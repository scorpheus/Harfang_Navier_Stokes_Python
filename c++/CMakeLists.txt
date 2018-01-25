cmake_minimum_required(VERSION 3.1)

project(navier_stokes)

set(CMAKE_CXX_STANDARD 11)

set(HARFANG_SDK "sdk" CACHE STRING "Harfang: Path to Harfang SDK")
get_filename_component(HARFANG_SDK ${HARFANG_SDK} REALPATH)

if(WIN32)
	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/Od /Ob0 /Zi")
else()
	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0")
endif()

if(CMAKE_CONFIGURATION_TYPES)
	set(HARFANG_LIBS_PATH ${HARFANG_SDK}/lib/${CMAKE_CFG_INTDIR})
	set(HARFANG_BIN_PATH ${HARFANG_SDK}/bin/${CMAKE_CFG_INTDIR})
else()
	set(HARFANG_LIBS_PATH ${HARFANG_SDK}/lib/${CMAKE_BUILD_TYPE})
	set(HARFANG_BIN_PATH ${HARFANG_SDK}/bin/${CMAKE_BUILD_TYPE})

	if(NOT EXISTS ${HARFANG_LIBS_PATH})
		message(FATAL_ERROR "Could not find Harfang libraries at " ${HARFANG_LIBS_PATH})
	endif()
endif()

message(STATUS "Harfang libraries at " ${HARFANG_LIBS_PATH})
message(STATUS "Harfang binaries at " ${HARFANG_BIN_PATH})

link_directories(${HARFANG_LIBS_PATH})
set(HARFANG_LIBS engine platform foundation lua)

if(WIN32)
	set(HARFANG_LIBS ${HARFANG_LIBS} Dbghelp Winmm Xinput dinput8 dxguid Shcore)
else()
	set(HARFANG_LIBS ${HARFANG_LIBS} X11 GL dl)
endif()

# PCH
function(configure_pch SRCS)
	if(CHIVAS_USE_PCH)
		set(PCH_SRC pch.cpp)
		set(PCH_HDR pch.h)
		if(MSVC)
			set_source_files_properties(${SRCS} PROPERTIES COMPILE_FLAGS "/Yu\"${PCH_HDR}\" /FI\"${PCH_HDR}\"")  
			set_source_files_properties(${PCH_SRC} PROPERTIES COMPILE_FLAGS "/Yc\"${PCH_HDR}\" /FI\"${PCH_HDR}\"")
		endif()
	endif()
endfunction()

if(CHIVAS_USE_PCH)
	if(MSVC)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Zm193")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Zm193")
	endif()
endif()

file(GLOB SRCS . *.cpp)
file(GLOB HDRS . *.h)

list(REMOVE_ITEM SRCS "${CMAKE_SOURCE_DIR}/navier_stokes_working_array2D.cpp")

add_executable(navier_stokes ${SRCS} ${HDRS})
configure_pch("${SRCS}")

target_include_directories(navier_stokes PUBLIC ${HARFANG_SDK}/include)
TARGET_LINK_LIBRARIES(navier_stokes ${HARFANG_LIBS})


add_custom_command(TARGET navier_stokes POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory ${HARFANG_SDK}/bin/${CMAKE_CFG_INTDIR} ${PROJECT_SOURCE_DIR}/executable
	COMMAND ${CMAKE_COMMAND} -E copy_directory ../assets ${PROJECT_SOURCE_DIR}/executable/assets
	COMMAND ${CMAKE_COMMAND} -E copy_directory $<TARGET_FILE_DIR:navier_stokes> "${PROJECT_SOURCE_DIR}/executable"
	COMMENT "Deploying navier_stokes")
	