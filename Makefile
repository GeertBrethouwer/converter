# make            build ./rsconv (any C++17 compiler; no libraries needed)
# make test       FFT test and round-trip tests (python3 with numpy)
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -pthread

rsconv: src/rsconv.cpp
	$(CXX) $(CXXFLAGS) src/rsconv.cpp -o $@

tests/test_fft: tests/test_fft.cpp src/rsconv.cpp
	$(CXX) $(CXXFLAGS) tests/test_fft.cpp -o $@

test: rsconv tests/test_fft
	./tests/test_fft
	python3 tests/test_roundtrip.py ./rsconv

clean:
	rm -f rsconv tests/test_fft

.PHONY: test clean
