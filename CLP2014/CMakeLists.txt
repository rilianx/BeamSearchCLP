# HelloWorld is a minimal sample creating, stepping and deleting a Bullet dynamics world

INCLUDE_DIRECTORIES(
../bullet-2.80-rev2531/src
)

LINK_DIRECTORIES(
../bullet-2.80-rev2531/src/BulletCollision
../bullet-2.80-rev2531/src/BulletDynamics
../bullet-2.80-rev2531/src/LinearMath
)

# add_subdirectory (../bullet-casa/bullet-2.79/src/BulletCollision)

LINK_LIBRARIES(
 BulletDynamics BulletCollision LinearMath 
)
    if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
    endif(COMMAND cmake_policy)


 SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# ADD_DEFINITION(-g -pg)

ADD_LIBRARY(_3dodp_lib Container.cpp Box.cpp Block.cpp CLP.cpp SortedList.cpp)

IF (WIN32)
	ADD_EXECUTABLE(MainCLP
		MainCLP-REVIEW.cpp 
		${BULLET_PHYSICS_SOURCE_DIR}/msvc/bullet.rc
	)
ELSE()

#~ ADD_EXECUTABLE(CLP-no-bpr
#~ MainCLP-no-bpr.cpp
#~ )
#~ target_link_libraries(CLP-no-bpr _3dodp_lib)
#~ 
#~ 
#~ ADD_EXECUTABLE(CLP-paper
#~ MainCLP-paper.cpp
#~ )
#~ target_link_libraries(CLP-paper _3dodp_lib)

ADD_EXECUTABLE(CLP2
MainCLP-REVIEW.cpp
)
target_link_libraries(CLP2 _3dodp_lib)

ENDIF()




IF (INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
			SET_TARGET_PROPERTIES(AppHelloWorld PROPERTIES  DEBUG_POSTFIX "_Debug")
			SET_TARGET_PROPERTIES(AppHelloWorld PROPERTIES  MINSIZEREL_POSTFIX "_MinsizeRel")
			SET_TARGET_PROPERTIES(AppHelloWorld PROPERTIES  RELWITHDEBINFO_POSTFIX "_RelWithDebugInfo")
ENDIF(INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
