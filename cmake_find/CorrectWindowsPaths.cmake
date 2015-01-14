# CorrectWindowsPaths - this module defines one macro
#
# CONVERT_CYGWIN_PATH( PATH )
#  This uses the command cygpath (provided by cygwin) to convert
#  unix-style paths into paths useable by cmake on windows

macro (CONVERT_CYGWIN_PATH _path)
  if (WIN32)
    if( NOT DEFINED CYGPATH_EXECUTABLE OR NOT EXISTS ${CYGPATH_EXECUTABLE} )
	  find_program(CYGPATH_EXECUTABLE NAMES cygpath HINTS
	            "/usr/bin/"
				"c:/cygwin/bin"
				"c:/cygwin64/bin"
				"${CYGWIN_INSTALL_PATH}/bin")
	  if(NOT EXISTS ${CYGPATH_EXECUTABLE}) 
        message(SEND_ERROR "cannot find cygpath.exe, please set CYGWIN_INSTALL_PATH variable with path where you have installed cygwin or CYGPATH_EXECUTABLE with path to cygpath.exe")
      endif()
	endif()
#	  message(STATUS "convert input")
#	  message(${${_path}})
	  string (REGEX REPLACE "\\?space\\?" " " ${_path} ${${_path}}) 
#	  message(STATUS "replace spaces")
#	  message(${${_path}})
      EXECUTE_PROCESS(COMMAND ${CYGPATH_EXECUTABLE} -m "${${_path}}"
        OUTPUT_VARIABLE ${_path})
      string (STRIP ${${_path}} ${_path})
	  string(REGEX REPLACE "\\\\" "/" ${_path} ${${_path}}) 
#	  message(STATUS "convert output")
#	  message(${${_path}})
	  
  endif (WIN32)
endmacro (CONVERT_CYGWIN_PATH)

