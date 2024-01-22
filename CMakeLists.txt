CMAKE_MINIMUM_REQUIRED(VERSION 3.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

INCLUDE(${CMAKE_SOURCE_DIR}/CCMacros.txt)  

# Definition of the project name
PROJECT(MultiStencilFastMarching)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
SET(Required_Packages
	 ITK
	 )

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

# Load boost program_options
#SET(BOOST_ROOT /home/export/personal/smercav/bin/boost_1_71_0/boost_1_71_0)
FIND_PACKAGE( Boost COMPONENTS program_options REQUIRED )



# Variables to be used in CMake
FOREACH(MI MultiStencilFastMarching MultiStencilFastMarching2D)
  SET (${MI}_SOURCE ${MI}.cxx)
  #SET (CMAKE_CXX_STANDARD 11)

  # add_executable -> creates a stand-alone executable to be run by console
  ADD_EXECUTABLE(${MI} ${${MI}_SOURCE})


  TARGET_LINK_LIBRARIES(${MI} 
	  ${ITK_LIBRARIES}
	  ${Boost_LIBRARIES}
	  ) 

ENDFOREACH(MI)

# Package the code
INSTALL(TARGETS MultiStencilFastMarching MultiStencilFastMarching2D RUNTIME DESTINATION bin)
INCLUDE(CPack)

