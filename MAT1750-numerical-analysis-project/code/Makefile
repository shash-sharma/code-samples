
MKFILE_PATH = ${CURDIR}

# Compiler
CXX = g++

OBJ_PATH = ./obj
INC_PATH = .
SRC_PATH = .
EXE_PATH = .
LIB_EXT = ${REPO_PATH}/External

COORD_DIR = ./common
MESH_DIR = ./mesh
MOM_DIR = ./MoM_integrals


# ============ Bessel functions ============

BESSEL_DIR = ./AmosBessel
BESSELobjects = ${BESSEL_DIR}/libAmosBessel.o ${BESSEL_DIR}/machcon.o ${BESSEL_DIR}/zbsubs.o


# ============ VTK ============

VTK_PATH = ${LIB_EXT}/VTK-8.1.1

VTK_LIB = -Wl,-rpath,${VTK_PATH}/build/lib -L${VTK_PATH}/build/lib
VTK_INC = -I${VTK_PATH}/Common/Core -I${VTK_PATH}/build/Common/Core -I${VTK_PATH}/Utilities/KWIML -I${VTK_PATH}/build/Utilities/KWIML -I${VTK_PATH}/Utilities/KWIML/vtkkwiml -I${VTK_PATH}/Common/DataModel -I${VTK_PATH}/build/Common/DataModel -I${VTK_PATH}/Common/ExecutionModel -I${VTK_PATH}/build/Common/ExecutionModel -I${VTK_PATH}/IO/XML -I${VTK_PATH}/build/IO/XML -I${VTK_PATH}/Common/Math -I${VTK_PATH}/build/Common/Math


# ============ Collections ============

ALLobjects = ${OBJ_PATH}/testMLFMA.o ${OBJ_PATH}/MLFMA.o ${OBJ_PATH}/MLFMA_tools.o ${OBJ_PATH}/structure.o ${OBJ_PATH}/VTK_visualization.o ${OBJ_PATH}/mesh.o ${BESSELobjects}


# ============ Compilation settings ============

CXXFLAGS = -O3 -std=c++11 -w -g -I${INC_PATH} -I${COORD_DIR} -I${MESH_DIR} -I${MOM_DIR}

LDFLAGS = -lgfortran -L${BESSEL_DIR} -lAmosBessel

LDLIBS = 
GPROF = 

COMPILE.cc = ${CXX} ${GPROF} ${CXXFLAGS}

SHARED_HEADERS = ${COORD_DIR}/coordinate_system.hpp

# Final executables
all: testMLFMA

# Main executable
testMLFMA: $(ALLobjects)
	${COMPILE.cc} ${LDFLAGS} -o testMLFMA $(ALLobjects) ${LDFLAGS}


${OBJ_PATH}/testMLFMA.o: ${SRC_PATH}/testMLFMA.cpp ${SHARED_HEADERS}
	${COMPILE.cc} -c ${SRC_PATH}/testMLFMA.cpp -o ${OBJ_PATH}/testMLFMA.o ${LDFLAGS}

${OBJ_PATH}/MLFMA.o: ${SRC_PATH}/MLFMA.cpp ${SHARED_HEADERS} ${INC_PATH}/MLFMA.hpp ${MESH_DIR}/mesh.hpp
	${COMPILE.cc} -c ${SRC_PATH}/MLFMA.cpp -o ${OBJ_PATH}/MLFMA.o ${LDFLAGS}

${OBJ_PATH}/MLFMA_tools.o: ${SRC_PATH}/MLFMA_tools.cpp ${SHARED_HEADERS} ${INC_PATH}/MLFMA_tools.hpp ${MESH_DIR}/mesh.hpp
	${COMPILE.cc} -c ${SRC_PATH}/MLFMA_tools.cpp -o ${OBJ_PATH}/MLFMA_tools.o ${LDFLAGS}

${OBJ_PATH}/structure.o: ${SRC_PATH}/structure.cpp ${SHARED_HEADERS} ${SRC_PATH}/structure.hpp ${MESH_DIR}/mesh.hpp
	${COMPILE.cc} -c ${SRC_PATH}/structure.cpp -o ${OBJ_PATH}/structure.o ${LDFLAGS}

${OBJ_PATH}/VTK_visualization.o: ${MESH_DIR}/VTK_visualization.cpp ${SHARED_HEADERS} ${MESH_DIR}/VTK_visualization.hpp ${MESH_DIR}/mesh.hpp
	${COMPILE.cc} -c ${MESH_DIR}/VTK_visualization.cpp -o ${OBJ_PATH}/VTK_visualization.o ${LDFLAGS}

${OBJ_PATH}/mesh.o: ${MESH_DIR}/mesh.cpp ${SHARED_HEADERS} ${MESH_DIR}/mesh.hpp
	${COMPILE.cc} -c ${MESH_DIR}/mesh.cpp -o ${OBJ_PATH}/mesh.o ${LDFLAGS}


# Clean
.PHONY: clean

clean:
	rm *.o -f main ${OBJ_PATH}/*.o

