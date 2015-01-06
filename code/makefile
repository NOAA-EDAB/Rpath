C_VER   = ecosim_r.cpp
H_VER   = ecosim_r.h
R_LIBS  = -lm 
R_OBJS  = ecosim_r.o
R_GPP   = g++
R_FLAGS = -shared -O3 -static-libgcc -static-libstdc++ -fPIC

all : ${R_OBJS}
	@echo "Making r verion"
	${R_GPP} -o ecosim.dll ${R_FLAGS} ${R_OBJS} ${R_LIBS}
  
ecosim_r.o : ${C_VER} ${H_VER}  
	${R_GPP} -c ${R_FLAGS} ${C_VER}
