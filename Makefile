# SLAMCore-8G top-level entry points.
#
#   make          configure + build
#   make test     run the CTest suite
#   make bench    regenerate the benchmark report under docs/
#   make clean    remove the build tree

BUILD_DIR ?= build
JOBS      ?= $(shell nproc 2>/dev/null || echo 4)
CMAKE     ?= cmake
BENCH_SEEDS ?= 1 2 3 4 5 6 7 8

.PHONY: all build test bench clean rebuild

all: build

$(BUILD_DIR)/CMakeCache.txt:
	$(CMAKE) -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release

build: $(BUILD_DIR)/CMakeCache.txt
	$(CMAKE) --build $(BUILD_DIR) -j$(JOBS)

test: build
	cd $(BUILD_DIR) && ctest --output-on-failure

# Regenerates every table in docs/benchmark_report.md from the built binary.
# Nothing in the report is hand-written, so the numbers can always be traced
# back to a command that reproduces them.
bench: build
	bash scripts/bench.sh $(BUILD_DIR)/slamcore $(BENCH_SEEDS)

clean:
	rm -rf $(BUILD_DIR)

rebuild: clean build
