#
# Copyright © Cambridge Numerical Solutions Ltd 2013
#
.SUFFIXES: .cu .cpp

all: ghost

TARGET = ghost

SRCDIR = src

SOURCE = ghost.cu

# Find picture and include dependencies from a tex file
%.d: %.cu
	nvcc $< --compiler-bindir=/usr/local -M | sed 's,\($*\)\.o[ :]*,\1 $@ : ,g' > $@

NVCC:= /usr/local/cuda-10.0/bin/nvcc
CUDA_DIR:=$(dir $(NVCC))../
CUDA_INC:=$(CUDA_DIR)/include
ARCH:=$(shell uname -m)

ifeq ($(ARCH),x86_64)
CUDA_LIB:=$(CUDA_DIR)/lib64
else
CUDA_LIB:=$(CUDA_DIR)/lib
endif

# ============================================================
# 从 ghost.hpp 提取网格参数, 编译前自动生成 eddie.cfg
# 用户只需修改 ghost.hpp 中的 N_CELLS_X/Y 和 L_DOMAIN_X/Y
# ============================================================
N_CELLS_X := $(shell grep 'N_CELLS_X =' ghost.hpp | sed 's/.*= \([0-9]*\).*/\1/')
N_CELLS_Y := $(shell grep 'N_CELLS_Y =' ghost.hpp | sed 's/.*= \([0-9]*\).*/\1/')
L_DOMAIN_X := $(shell grep 'L_DOMAIN_X =' ghost.hpp | sed 's/.*= \([0-9.]*\).*/\1/')
L_DOMAIN_Y := $(shell grep 'L_DOMAIN_Y =' ghost.hpp | sed 's/.*= \([0-9.]*\).*/\1/')

eddie.cfg: eddie.cfg.template ghost.hpp
	@echo "# Generating eddie.cfg: N_CELLS=($(N_CELLS_X),$(N_CELLS_Y)) L_DOMAIN=($(L_DOMAIN_X),$(L_DOMAIN_Y))"
	sed -e 's/@N_CELLS_X@/$(N_CELLS_X)/g' \
	    -e 's/@N_CELLS_Y@/$(N_CELLS_Y)/g' \
	    -e 's/@L_DOMAIN_X@/$(L_DOMAIN_X)/g' \
	    -e 's/@L_DOMAIN_Y@/$(L_DOMAIN_Y)/g' \
	    $< > $@

$(TARGET): $(SOURCE) eddie.cfg
	$(NVCC) $< --expt-relaxed-constexpr -o $@ --compiler-bindir=/usr/bin/ -Xptxas -v -arch=compute_60 -code=sm_60 -L/usr/local/cuda-10.0/lib64 -L/usr/local/cuda-10.0/samples/common/lib/linux/ -I/usr/local/cuda-10.0/samples/common/inc -I../../src -L$(HOME)/boost_1_58_0/stage/lib -I$(HOME)/boost_1_58_0 -I. -lrt -lcuda -lboost_thread -lboost_system -lglut -lpng -lconfig++
	
	#$(NVCC) $< -o $@ --compiler-bindir=/usr/bin/ -Xptxas -v -arch=compute_60 -code=sm_60 -L/usr/local/cuda-10.0/lib64 -L/usr/local/cuda-10.0/samples/common/lib/linux/ -I/usr/local/cuda-10.0/samples/common/inc -I../../src -I. -lrt -lcuda -lboost_thread -lboost_system -lglut -lpng -lconfig++
    #$(NVCC) $< -o $@ -std=c++11 -Xptxas -v -arch=compute_60 -code=sm_60 -L$(CUDA_LIB) -L$(CUDA_DIR)/samples/common/lib/linux/ -I$(CUDA_DIR)/samples/common/inc -I../../src -L$(HOME)/local/boost_1_58_0/stage/lib -I$(HOME)/local/boost_1_58_0 -I. -lrt -lboost_thread -lboost_system -lcuda -lglut -lpng -lconfig++