# Search program.  This is NOT part of the proof; see README.txt.
CAPD_DIR ?= $(HOME)/.local/capd
CAPD_CONFIG = $(CAPD_DIR)/bin/capd-config

CXXFLAGS = $(shell $(CAPD_CONFIG) --cflags)
LIBS     = $(shell $(CAPD_CONFIG) --libs)

all: certify_range

certify_range: certify_range.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(LIBS)

clean:
	rm -f certify_range

.PHONY: all clean
