# make            build ./rsconv and ./rsspec (any C++17 compiler; no libraries needed)
# make test       FFT test, round-trip tests and spectra tests (python3 with numpy)
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -pthread

all: rsconv rsspec

rsconv: src/rsconv.cpp src/rsio.hpp
	$(CXX) $(CXXFLAGS) src/rsconv.cpp -o $@

rsspec: src/rsspec.cpp src/rsio.hpp
	$(CXX) $(CXXFLAGS) src/rsspec.cpp -o $@

tests/test_fft: tests/test_fft.cpp src/rsconv.cpp src/rsio.hpp
	$(CXX) $(CXXFLAGS) tests/test_fft.cpp -o $@

test: rsconv rsspec tests/test_fft
	./tests/test_fft
	python3 tests/test_roundtrip.py ./rsconv
	python3 tests/test_spectra.py .

clean:
	rm -f rsconv rsspec tests/test_fft

.PHONY: all test clean
