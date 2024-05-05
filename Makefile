# Makefile for SCPD_Project

# Variables
PROJECT_BUILD_DIR := cmake-build
TARGET := SCPD_Project
MAKE_PROGRAM := $(shell which ninja)
BENCHMARK_BIN_DIR := benchmark/bin

# Default target
all: release debug genJson suiteTest

# Release target
release:
	$(MAKE) build BUILD_TYPE=Release

# Debug target
debug:
	$(MAKE) build BUILD_TYPE=Debug

# Build rule
build:
	mkdir -p $(PROJECT_BUILD_DIR)/$(BUILD_TYPE) && \
	cmake -S . -B $(PROJECT_BUILD_DIR)/$(BUILD_TYPE) \
		-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
		-DCMAKE_MAKE_PROGRAM=$(MAKE_PROGRAM) \
		-G Ninja && \
	cmake --build $(PROJECT_BUILD_DIR)/$(BUILD_TYPE) --target $(TARGET) -j 10

# Generate JSON rule
genJson: $(BENCHMARK_BIN_DIR)/genJson
$(BENCHMARK_BIN_DIR)/genJson: benchmark/src/genJson.cc
	mkdir -p $(BENCHMARK_BIN_DIR)
	g++ -o $@ $<

# Suite test rule
suiteTest: $(BENCHMARK_BIN_DIR)/suiteTest
$(BENCHMARK_BIN_DIR)/suiteTest: benchmark/src/suiteTest.cc
	mkdir -p $(BENCHMARK_BIN_DIR)
	g++ -o $@ $<

sequentialAP: $(PROJECT_BUILD_DIR)/seuentialAP
$(PROJECT_BUILD_DIR)/suiteTest: src/data_structures/sequential_ap.cc src/main.cc include/utilities/nbody_helpers.h include/data_structures/sequential_ap.h
	mkdir -p $(PROJECT_BUILD_DIR)
	g++ -o $@ $<
# Clean rule
clean:
	rm -rf $(PROJECT_BUILD_DIR) $(BENCHMARK_BIN_DIR)

.PHONY: all release debug genJson suiteTest build clean