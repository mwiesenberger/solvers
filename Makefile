device=omp
FELTOR_PATH=../feltor

#configure machine
include $(FELTOR_PATH)/config/default.mk
include $(FELTOR_PATH)/config/version.mk
include $(FELTOR_PATH)/config/*.mk
include $(FELTOR_PATH)/config/devices/devices.mk

INCLUDE+=-I$(FELTOR_PATH)/inc/

all: solvers

solvers: solvers.cpp
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -g

.PHONY: clean

clean:
	rm -f solvers
