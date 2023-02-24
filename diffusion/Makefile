device=gpu
FELTOR_PATH=../../feltor

#configure machine
include $(FELTOR_PATH)/config/default.mk
include $(FELTOR_PATH)/config/version.mk
include $(FELTOR_PATH)/config/*.mk
include $(FELTOR_PATH)/config/devices/devices.mk

INCLUDE+=-I$(FELTOR_PATH)/inc/


INCLUDE+= -I../../inc   # other project libraries

all: diffusion

diffusion: diffusion.cpp diffusion.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) $(LAPACKLIB) -g

.PHONY: clean

clean:
	rm -rf diffusion
