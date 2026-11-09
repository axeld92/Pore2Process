CXX := /opt/homebrew/bin/g++-15
CC := /opt/homebrew/bin/gcc-15
CXXFLAGS := -O2 -Wall -std=c++20 -I/opt/homebrew/include
LDFLAGS := -L$(CURDIR)/core -L/opt/homebrew/lib
PARALLEL_MODE := OFF
PLATFORMS := CPU_SISD
FLOATING_POINT_TYPE := double
FEATURES :=
USE_EMBEDDED_DEPENDENCIES := OFF
