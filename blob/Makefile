device=gpu
FELTOR_PATH=../../feltor

#configure machine
include $(FELTOR_PATH)/config/default.mk
include $(FELTOR_PATH)/config/version.mk
include $(FELTOR_PATH)/config/*.mk
include $(FELTOR_PATH)/config/devices/devices.mk

INCLUDE+=-I$(FELTOR_PATH)/inc/


INCLUDE+= -I../../inc   # other project libraries

all: toefl toefl_hpc toefl_mpi

toefl: toefl.cpp toefl.h parameters.h diag.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) -DWITH_GLFW -g

toefl_hpc: toefl.cpp toefl.h parameters.h diag.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) -g

toefl_mpi: toefl.cpp toefl.h parameters.h diag.h
	$(MPICC) $(OPT) $(MPICFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) $(VERSION_FLAGS)

doc:
	pdflatex -shell-escape ./toefl.tex;
	bibtex toefl.aux;
	pdflatex -shell-escape ./toefl.tex;
	pdflatex -shell-escape ./toefl.tex;


.PHONY: clean doc

clean:
	rm -rf toefl toefl_hpc toefl_mpi toefl.aux toefl.log toefl.out toefl.pyg toefl.pdf toefl.bbl toefl.blg
