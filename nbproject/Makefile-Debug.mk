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
FC=f95
AS=as

# Macros
CND_PLATFORM=Cygwin-Windows
CND_DLIB_EXT=dll
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/FIBR3D_app.o \
	${OBJECTDIR}/FillingAlgorithms.o \
	${OBJECTDIR}/Inspection.o \
	${OBJECTDIR}/MMDS.o \
	${OBJECTDIR}/PartsOrientation.o \
	${OBJECTDIR}/admesh/connect.o \
	${OBJECTDIR}/admesh/normals.o \
	${OBJECTDIR}/admesh/shared.o \
	${OBJECTDIR}/admesh/stl_io.o \
	${OBJECTDIR}/admesh/stlinit.o \
	${OBJECTDIR}/admesh/util.o \
	${OBJECTDIR}/clipper.o \
	${OBJECTDIR}/clipper_utils.o \
	${OBJECTDIR}/edges.o \
	${OBJECTDIR}/export.o \
	${OBJECTDIR}/lhs/geneticLHS.o \
	${OBJECTDIR}/lhs/improvedLHS.o \
	${OBJECTDIR}/lhs/maximinLHS.o \
	${OBJECTDIR}/lhs/optSeededLHS.o \
	${OBJECTDIR}/lhs/optimumLHS.o \
	${OBJECTDIR}/lhs/randomLHS.o \
	${OBJECTDIR}/lhs/utilityLHS.o \
	${OBJECTDIR}/my_stl.o \
	${OBJECTDIR}/optimization.o \
	${OBJECTDIR}/pattern.o \
	${OBJECTDIR}/pswarm.o \
	${OBJECTDIR}/splines.o \
	${OBJECTDIR}/triangleIntersection/NoDivTriTriIsect.o


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
LDLIBSOPTIONS=-Lglpk lp_solve/lpsolve55.dll blas/tmglib_LINUX.a blas/lapack_LINUX.a blas/blas_LINUX.a -lf2c

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fibr3d_app.exe
	${CP} lp_solve/lpsolve55.dll ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fibr3d_app.exe: lp_solve/lpsolve55.dll

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fibr3d_app.exe: blas/tmglib_LINUX.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fibr3d_app.exe: blas/lapack_LINUX.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fibr3d_app.exe: blas/blas_LINUX.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fibr3d_app.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fibr3d_app ${OBJECTFILES} ${LDLIBSOPTIONS} -Wl,--enable-auto-import

${OBJECTDIR}/FIBR3D_app.o: FIBR3D_app.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FIBR3D_app.o FIBR3D_app.cpp

${OBJECTDIR}/FillingAlgorithms.o: FillingAlgorithms.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FillingAlgorithms.o FillingAlgorithms.cpp

${OBJECTDIR}/Inspection.o: Inspection.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Inspection.o Inspection.cpp

${OBJECTDIR}/MMDS.o: MMDS.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MMDS.o MMDS.cpp

${OBJECTDIR}/PartsOrientation.o: PartsOrientation.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/PartsOrientation.o PartsOrientation.cpp

${OBJECTDIR}/admesh/connect.o: admesh/connect.c
	${MKDIR} -p ${OBJECTDIR}/admesh
	${RM} "$@.d"
	$(COMPILE.c) -g -DHAVE_CONFIG_H=1 -IArmadillo -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/admesh/connect.o admesh/connect.c

${OBJECTDIR}/admesh/normals.o: admesh/normals.c
	${MKDIR} -p ${OBJECTDIR}/admesh
	${RM} "$@.d"
	$(COMPILE.c) -g -DHAVE_CONFIG_H=1 -IArmadillo -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/admesh/normals.o admesh/normals.c

${OBJECTDIR}/admesh/shared.o: admesh/shared.c
	${MKDIR} -p ${OBJECTDIR}/admesh
	${RM} "$@.d"
	$(COMPILE.c) -g -DHAVE_CONFIG_H=1 -IArmadillo -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/admesh/shared.o admesh/shared.c

${OBJECTDIR}/admesh/stl_io.o: admesh/stl_io.c
	${MKDIR} -p ${OBJECTDIR}/admesh
	${RM} "$@.d"
	$(COMPILE.c) -g -DHAVE_CONFIG_H=1 -IArmadillo -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/admesh/stl_io.o admesh/stl_io.c

${OBJECTDIR}/admesh/stlinit.o: admesh/stlinit.c
	${MKDIR} -p ${OBJECTDIR}/admesh
	${RM} "$@.d"
	$(COMPILE.c) -g -DHAVE_CONFIG_H=1 -IArmadillo -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/admesh/stlinit.o admesh/stlinit.c

${OBJECTDIR}/admesh/util.o: admesh/util.c
	${MKDIR} -p ${OBJECTDIR}/admesh
	${RM} "$@.d"
	$(COMPILE.c) -g -DHAVE_CONFIG_H=1 -IArmadillo -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/admesh/util.o admesh/util.c

${OBJECTDIR}/clipper.o: clipper.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/clipper.o clipper.cpp

${OBJECTDIR}/clipper_utils.o: clipper_utils.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/clipper_utils.o clipper_utils.cpp

${OBJECTDIR}/edges.o: edges.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/edges.o edges.cpp

${OBJECTDIR}/export.o: export.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/export.o export.cpp

${OBJECTDIR}/lhs/geneticLHS.o: lhs/geneticLHS.cpp
	${MKDIR} -p ${OBJECTDIR}/lhs
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lhs/geneticLHS.o lhs/geneticLHS.cpp

${OBJECTDIR}/lhs/improvedLHS.o: lhs/improvedLHS.cpp
	${MKDIR} -p ${OBJECTDIR}/lhs
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lhs/improvedLHS.o lhs/improvedLHS.cpp

${OBJECTDIR}/lhs/maximinLHS.o: lhs/maximinLHS.cpp
	${MKDIR} -p ${OBJECTDIR}/lhs
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lhs/maximinLHS.o lhs/maximinLHS.cpp

${OBJECTDIR}/lhs/optSeededLHS.o: lhs/optSeededLHS.cpp
	${MKDIR} -p ${OBJECTDIR}/lhs
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lhs/optSeededLHS.o lhs/optSeededLHS.cpp

${OBJECTDIR}/lhs/optimumLHS.o: lhs/optimumLHS.cpp
	${MKDIR} -p ${OBJECTDIR}/lhs
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lhs/optimumLHS.o lhs/optimumLHS.cpp

${OBJECTDIR}/lhs/randomLHS.o: lhs/randomLHS.cpp
	${MKDIR} -p ${OBJECTDIR}/lhs
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lhs/randomLHS.o lhs/randomLHS.cpp

${OBJECTDIR}/lhs/utilityLHS.o: lhs/utilityLHS.cpp
	${MKDIR} -p ${OBJECTDIR}/lhs
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/lhs/utilityLHS.o lhs/utilityLHS.cpp

${OBJECTDIR}/my_stl.o: my_stl.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/my_stl.o my_stl.cpp

${OBJECTDIR}/optimization.o: optimization.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/optimization.o optimization.cpp

${OBJECTDIR}/pattern.o: pattern.c
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -DHAVE_CONFIG_H=1 -IArmadillo -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/pattern.o pattern.c

${OBJECTDIR}/pswarm.o: pswarm.c
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -DHAVE_CONFIG_H=1 -IArmadillo -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/pswarm.o pswarm.c

${OBJECTDIR}/splines.o: splines.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Ilp_solve -IArmadillo -Ilhs -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/splines.o splines.cpp

${OBJECTDIR}/triangleIntersection/NoDivTriTriIsect.o: triangleIntersection/NoDivTriTriIsect.c
	${MKDIR} -p ${OBJECTDIR}/triangleIntersection
	${RM} "$@.d"
	$(COMPILE.c) -g -DHAVE_CONFIG_H=1 -IArmadillo -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/triangleIntersection/NoDivTriTriIsect.o triangleIntersection/NoDivTriTriIsect.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} -r ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/lpsolve55.dll
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fibr3d_app.exe

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
