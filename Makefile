CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
SDL_CFLAGS := $(shell sdl2-config --cflags 2>/dev/null)
SDL_LIBS   := $(shell sdl2-config --libs 2>/dev/null)

all: double_pendulum analysis test_physics

double_pendulum: main.cpp pendulum.hpp
	$(CXX) $(CXXFLAGS) $(SDL_CFLAGS) main.cpp -o $@ $(SDL_LIBS) -lm

analysis: analysis.cpp pendulum.hpp
	$(CXX) $(CXXFLAGS) analysis.cpp -o $@

test_physics: test_physics.cpp pendulum.hpp
	$(CXX) $(CXXFLAGS) test_physics.cpp -o $@

test: test_physics
	./test_physics

figures: analysis
	./analysis figures

clean:
	rm -f double_pendulum analysis test_physics capture.bmp

.PHONY: all test figures clean
