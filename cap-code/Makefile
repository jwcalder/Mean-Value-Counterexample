# Requires CAPD (https://github.com/CAPDGroup/CAPD) built and installed.
# Set CAPD_DIR to the install prefix if it is not ~/.local/capd.
CAPD_DIR ?= $(HOME)/.local/capd
CAPD_CONFIG = $(CAPD_DIR)/bin/capd-config

CXXFLAGS = $(shell $(CAPD_CONFIG) --cflags)
LIBS     = $(shell $(CAPD_CONFIG) --libs)

all: verify

verify: verify.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(LIBS)

clean:
	rm -f verify

.PHONY: all clean
