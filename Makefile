#---------------------------------
BIN =
SRC_TEST = test.f90
LIB = libgreat.a
SRC1 = module_mode_integrator.f90 module_mode_analysis.f90
SRC2 = module_background.f90 module_eigen.f90 
SRC3 = 	module_param.f90 module_fast_inv.f90
DEPEND = 
#-----------------------------------


ROOT_PATH =  $(shell /bin/pwd)
MODULES_PATH = $(shell /bin/pwd)
include $(ROOT_PATH)/local_settings
ifneq ("$(wildcard  $(shell /bin/pwd)/very_local_settings)","")
 include $(shell /bin/pwd)/very_local_settings
endif
include $(ROOT_PATH)/makefile.inc

