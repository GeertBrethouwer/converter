# make            build ./rsconv and ./rsspec (any C++17 compiler; no libraries needed)
# make test       FFT test, round-trip tests and spectra tests (python3 with numpy)
#
# The binaries run only on the processor type they were built for (e.g. x86_64 on a login node, aarch64 on the
# GH200 nodes of Arrhenius): make rebuilds them when it is run on another type (build/.arch). The job scripts run
# make on the node they run on, so they always have binaries for it.
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -pthread
ARCH     := $(shell uname -m)
STAMP    := build/.arch

all: rsconv rsspec

$(STAMP): FORCE
	@mkdir -p build
	@echo '$(ARCH) $(CXX) $(CXXFLAGS)' | cmp -s - $@ || echo '$(ARCH) $(CXX) $(CXXFLAGS)' > $@

rsconv: src/rsconv.cpp src/rsio.hpp $(STAMP)
	$(CXX) $(CXXFLAGS) src/rsconv.cpp -o $@

rsspec: src/rsspec.cpp src/rsio.hpp $(STAMP)
	$(CXX) $(CXXFLAGS) src/rsspec.cpp -o $@

tests/test_fft: tests/test_fft.cpp src/rsconv.cpp src/rsio.hpp $(STAMP)
	$(CXX) $(CXXFLAGS) tests/test_fft.cpp -o $@

test: rsconv rsspec tests/test_fft
	./tests/test_fft
	python3 tests/test_roundtrip.py ./rsconv
	python3 tests/test_spectra.py .

clean:
	rm -rf rsconv rsspec tests/test_fft build

FORCE:

.PHONY: all test clean FORCE
