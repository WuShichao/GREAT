#************************************************************************
# Makefile include
#************************************************************************
#
# GPL  License:
# -------------
#
#    GREAT = General Relativistic Eigenmode Analysis Tool
#    Copyright (C) 2018 Alejandro Torres-Forne and Pablo Cerda-Duran
#
#    This file is part of GREAT.
#
#    GREAT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    GREAT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GREAT.  If not, see <http://www.gnu.org/licenses/>.
#
#************************************************************************

MAIN_PATH = $(shell /bin/pwd)

SRC_PATH = $(MAIN_PATH)/src
OBJ_PATH = $(MAIN_PATH)/obj
MOD_PATH = $(MAIN_PATH)/include
INC_PATH = $(MAIN_PATH)/include
LIB_PATH = $(MAIN_PATH)/lib
BIN_PATH = $(MAIN_PATH)/bin
DOC_PATH = $(MAIN_PATH)/doc

ifdef BIN
 EXEC = $(BIN_PATH)/$(BIN)
else
 EXEC =
endif

ifndef FCPP_FLAG
FCPP_FLAG = -cpp
endif


vpath %.f90 $(SRC_PATH)
vpath %.F $(SRC_PATH)
vpath %.C $(SRC_PATH)
vpath %.o $(OBJ_PATH)
vpath %.mod $(MOD_PATH)
vpath %.a $(LIB_PATH)
vpath %.tag $(TAG_PATH)

OBJ_TEST_F90 := $(filter %.o, $(SRC_TEST:.f90=.o))
OBJ1_F90 := $(filter %.o, $(SRC1:.f90=.o))
OBJ2_F90 := $(filter %.o, $(SRC2:.f90=.o))
OBJ3_F90 := $(filter %.o, $(SRC3:.f90=.o))
OBJ_TEST_F77 := $(filter %.o, $(SRC_TEST:.F=.o))
OBJ1_F77 := $(filter %.o, $(SRC1:.F=.o))
OBJ2_F77 := $(filter %.o, $(SRC2:.F=.o))
OBJ3_F77 := $(filter %.o, $(SRC3:.F=.o))
OBJ_TEST_CC := $(filter %.o, $(SRC_TEST:.C=.o))
OBJ1_CC := $(filter %.o, $(SRC1:.C=.o))
OBJ2_CC := $(filter %.o, $(SRC2:.C=.o))
OBJ3_CC := $(filter %.o, $(SRC3:.C=.o))

OBJ_TEST = $(OBJ_TEST_F90) $(OBJ_TEST_F77) $(OBJ_TEST_CC)
ifdef SRC_TEST
 EXEC_TEST_F90 := $(patsubst %.f90,%, $(filter %.f90, $(SRC_TEST)))
 EXEC_TEST_F77 := $(patsubst %.F,%, $(filter %.F, $(SRC_TEST)))
 EXEC_TEST_CC := $(patsubst %.C,%, $(filter %.C, $(SRC_TEST)))
 EXEC_TEST := $(EXEC_TEST_F90) $(EXEC_TEST_F77) $(EXEC_TEST_CC) 
else
 EXEC_TEST =
endif

SRC := $(SRC_TEST) $(SRC1) $(SRC2) $(SRC3)	
OBJ_F90 := $(OBJ1_F90) $(OBJ2_F90) $(OBJ3_F90)
OBJ_F77 := $(OBJ1_F77) $(OBJ2_F77) $(OBJ3_F77)
OBJ_CC := $(OBJ1_CC) $(OBJ2_CC) $(OBJ3_CC) $(OBJ)

#double substitution
expand = $(MODULES_PATH)/$(1)/lib/lib$(1).a
DEPEND_LIB = $(foreach x, $(DEPEND), $(call expand,$(x)))

INC_PATHS = $(INC_PATH)$(DEPEND:%= -I$(MODULES_PATH)/%/include)

all: $(LIB) $(EXEC)

test: $(LIB) $(EXEC) $(EXEC_TEST)


$(EXEC): $(DEPEND_LIB) $(OBJ_F90) $(OBJ_F77) $(OBJ_CC) $(LIB) 
	@echo "---------------------- Building executable:" $@ "-------------------------"
	cd $(OBJ_PATH); $(LD) $(LD_FLAGS) -o $(EXEC) $(OBJ_F90) $(OBJ_F77) $(OBJ_CC) $(DEPEND_LIB) $(CC_LIB) $(FC_LIB) $(LD_LIB) -I$(INC_PATH) 

$(EXEC_TEST):%:%.o $(DEPEND_LIB) $(OBJ_F90) $(OBJ_F77) $(OBJ_CC) $(LIB) 
	@echo "---------------------- Building executable (test):" $@ "-------------------------"
	cd $(OBJ_PATH); $(LD) $(LD_FLAGS) -o $(BIN_PATH)/$@ $(OBJ_PATH)/$@.o $(OBJ_F90) $(OBJ_F77) $(OBJ_CC) $(DEPEND_LIB) $(CC_LIB) $(FC_LIB) $(LD_LIB) -I$(INC_PATH) 

$(LIB):%.a: $(DEPEND_LIB) $(OBJ_F90) $(OBJ_F77) $(OBJ_CC)
	@echo "---------------------- Building library:" $@ "-------------------------"
	cd $(OBJ_PATH); $(AR) $(LIB_PATH)/$@ $(OBJ_F90) $(OBJ_F77) $(OBJ_CC)
#	cd $(LIB_PATH); $(RANLIB) $@ 

$(OBJ_TEST_F90):%.o:%.f90 $(DEPEND_LIB) $(OBJ1_F90) $(OBJ2_F90) $(OBJ3_F90)
	@echo "---------------------- Compiling level test code objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(FC) $(FCPP_FLAG) -c $(FC_FLAGS) $(FC_INCLUDES) -I$(INC_PATHS) $(MOD_FLAG)$(MOD_PATH)  $<

$(OBJ1_F90):%.o:%.f90 $(DEPEND_LIB) $(OBJ2_F90) $(OBJ3_F90)
	@echo "---------------------- Compiling level 1 objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(FC) $(FCPP_FLAG) -c $(FC_FLAGS) $(FC_INCLUDES) -I$(INC_PATHS) $(MOD_FLAG)$(MOD_PATH)  $<

$(OBJ2_F90):%.o:%.f90 $(DEPEND_LIB) $(OBJ3_F90)
	@echo "---------------------- Compiling level 2 objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(FC) $(FCPP_FLAG) -c $(FC_FLAGS) $(FC_INCLUDES) -I$(INC_PATHS) $(MOD_FLAG)$(MOD_PATH)  $<

$(OBJ3_F90):%.o:%.f90 $(DEPEND_LIB)
	@echo "---------------------- Compiling level 3 objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(FC) $(FCPP_FLAG) -c $(FC_FLAGS) $(FC_INCLUDES) -I$(INC_PATHS) $(MOD_FLAG)$(MOD_PATH)  $<

$(OBJ_TEST_F77):%.o:%.F $(DEPEND_LIB) $(OBJ1_F77) $(OBJ2_F77) $(OBJ3_F77)
	@echo "---------------------- Compiling level test code objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(FC) $(FCPP_FLAG) -c $(FC_FLAGS) $(FC_INCLUDES) -I$(INC_PATHS) $(MOD_FLAG)$(MOD_PATH)  $<

$(OBJ1_F77):%.o:%.F $(DEPEND_LIB) $(OBJ2_F77) $(OBJ3_F77)
	@echo "---------------------- Compiling level 1 objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(FC) $(FCPP_FLAG) -c $(FC_FLAGS) $(FC_INCLUDES) -I$(INC_PATHS) $(MOD_FLAG)$(MOD_PATH)  $<

$(OBJ2_F77):%.o:%.F $(DEPEND_LIB) $(OBJ3_F77)
	@echo "---------------------- Compiling level 2 objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(FC) $(FCPP_FLAG) -c $(FC_FLAGS) $(FC_INCLUDES) -I$(INC_PATHS) $(MOD_FLAG)$(MOD_PATH)  $<

$(OBJ3_F77):%.o:%.F $(DEPEND_LIB)
	@echo "---------------------- Compiling level 3 objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(FC) $(FCPP_FLAG) -c $(FC_FLAGS) $(FC_INCLUDES) -I$(INC_PATHS) $(MOD_FLAG)$(MOD_PATH)  $<

$(OBJ_TEST_CC):%.o:%.C $(DEPEND_LIB) $(OBJ1_CC) $(OBJ2_CC) $(OBJ3_CC)
	@echo "---------------------- Compiling level 1 objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(CC)  -c $(CC_FLAGS) $(CC_INCLUDES)  $<

$(OBJ1_CC):%.o:%.C $(DEPEND_LIB) $(OBJ2_CC) $(OBJ3_CC)
	@echo "---------------------- Compiling level 1 objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(CC)  -c $(CC_FLAGS) $(CC_INCLUDES)  $<

$(OBJ2_CC):%.o:%.C $(DEPEND_LIB) $(OBJ3_CC)
	@echo "---------------------- Compiling level 2 objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(CC)  -c $(CC_FLAGS) $(CC_INCLUDES)  $<
	     
$(OBJ3_CC):%.o:%.C $(DEPEND_LIB)
	@echo "---------------------- Compiling level 3 objects:" $@ "-------------------------"
	cd $(OBJ_PATH); $(CC)  -c $(CC_FLAGS) $(CC_INCLUDES)  $<


.PHONY: debug always

always:

debug:
	@echo "MAIN_PATH   : "$(MAIN_PATH)
	@echo "SRC_PATH    : "$(SRC_PATH)
	@echo "OBJ_PATH    : "$(OBJ_PATH)
	@echo "DEPEND      : "$(DEPEND)
	@echo "DEPEND_LIB  : "$(DEPEND_LIB)
	@echo "INC_PATHS   : "$(INC_PATHS)

doc:$(DOC_PATH)/doc.tag 

$(DOC_PATH)/doc.tag:$(SRC)
	cat ../../doc/doxyconf.coco doc/doxyconf.local | doxygen -  

clean:
	cd $(OBJ_PATH); /bin/rm -f *.o
	cd $(MOD_PATH); /bin/rm -f *.mod
	cd $(LIB_PATH); /bin/rm -f *.a
	cd $(BIN_PATH); /bin/rm -f $(BIN) $(EXEC_TEST)
#	cd $(DOC_PATH); /bin/rm -f *.tag ; /bin/rm -fr html
