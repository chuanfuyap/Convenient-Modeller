#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/SSsolverbridge.o \
	${OBJECTDIR}/TCsolverbridge.o \
	${OBJECTDIR}/TCsolverbridgePK.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/usr/local/lib/ -lODES -lsbml

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/SOSlibJNIC.so

dist/SOSlibJNIC.so: ${OBJECTFILES}
	${MKDIR} -p dist
	${LINK.c} -o dist/SOSlibJNIC.so ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/SSsolverbridge.o: SSsolverbridge.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SSsolverbridge.o SSsolverbridge.c

${OBJECTDIR}/TCsolverbridge.o: TCsolverbridge.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/TCsolverbridge.o TCsolverbridge.c

${OBJECTDIR}/TCsolverbridgePK.o: TCsolverbridgePK.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/TCsolverbridgePK.o TCsolverbridgePK.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/SOSlibJNIC.so

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
