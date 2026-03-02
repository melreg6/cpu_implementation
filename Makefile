# ================================
#   Makefile for LJ MD Engine
#   With Pthreads Support (CPU only)
# ================================

CC      = gcc
CFLAGS  = -O2 -std=c11 -Wall -Wextra -Iinclude -pthread
LDFLAGS = -lm -lpthread

CXX     = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra

PYTHON  = /opt/homebrew/bin/python3     # for plot_lj.py

# Number of threads for pthread neighbor list (default: 4)
# Override with: make NTHREADS=8
NTHREADS ?= 4
CFLAGS += -DNBL_NUM_THREADS=$(NTHREADS)

SRC = \
    src/main.c \
    src/md.c \
    src/cell_list.c \
    src/neighbor_list.c \
    src/neighbor_list_pthread.c \
    src/md_io.c \
    src/pdb_importer.c \
    src/cell_mt.c \
    src/cell_manhattan_pthread.c

OBJ = $(SRC:.c=.o)

TARGET      = md_run
COMPARE_BIN = compare_pdb
COMPARE_SRC = compare.cpp

# -----------------------------------------
# Default build
# -----------------------------------------
all: output $(TARGET)
	@echo ""
	@echo "Build complete with $(NTHREADS) threads for pthread neighbor list"
	@echo "To change: make clean && make NTHREADS=8"
	@echo ""

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET) $(LDFLAGS)

# -----------------------------------------
# Build comparison binary
# -----------------------------------------
$(COMPARE_BIN): $(COMPARE_SRC)
	$(CXX) $(CXXFLAGS) $(COMPARE_SRC) -o $(COMPARE_BIN)

# -----------------------------------------
# Create output directory automatically
# -----------------------------------------
output:
	mkdir -p output

# -----------------------------------------
# Build object files
# -----------------------------------------
src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

# -----------------------------------------
# Run simulation with default PDB
# -----------------------------------------
run: $(TARGET) output
	./md_run input.pdb

# -----------------------------------------
# Run comparison tests (writes to text file)
# -----------------------------------------
compare: $(COMPARE_BIN)
	@echo "Writing comparison report to output/compare_report.txt"
	@rm -f output/compare_report.txt

	@echo "========================================" >> output/compare_report.txt
	@echo "PDB COMPARISON REPORT" >> output/compare_report.txt
	@echo "Energy Conservation & Numerical Stability Test" >> output/compare_report.txt
	@echo "========================================" >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "1) Compare all methods against each other (after autotuning steps)" >> output/compare_report.txt
	@echo "----------------------------------------" >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL vs CELL-LIST:" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/cell_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL vs NEIGHBOR-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL vs NEIGHBOR-LIST (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST vs NEIGHBOR-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_positions.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST vs NEIGHBOR-LIST (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "NEIGHBOR-LIST (Serial) vs NEIGHBOR-LIST (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/nbl_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "2) Final positions vs all autotuning outputs" >> output/compare_report.txt
	@echo "----------------------------------------" >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs FULL:" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/full_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs CELL-LIST:" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/cell_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs NEIGHBOR-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs NEIGHBOR-LIST (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "3) Input vs all outputs (drift from initial)" >> output/compare_report.txt
	@echo "----------------------------------------" >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs FULL:" >> output/compare_report.txt
	@./$(COMPARE_BIN) input.pdb output/full_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs CELL-LIST:" >> output/compare_report.txt
	@./$(COMPARE_BIN) input.pdb output/cell_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs NEIGHBOR-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) input.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs NEIGHBOR-LIST (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) input.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs FINAL:" >> output/compare_report.txt
	@./$(COMPARE_BIN) input.pdb output/final_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "========================================" >> output/compare_report.txt
	@echo "Report complete." >> output/compare_report.txt
	@echo "Comparison report written to output/compare_report.txt"
	@cat output/compare_report.txt

# -----------------------------------------
# Plot energy graphs (Python script)
# -----------------------------------------
plot:
	$(PYTHON) test_scripts/plot_lj.py

# -----------------------------------------
# Basic clean
# -----------------------------------------
clean:
	rm -f $(OBJ) $(TARGET) $(COMPARE_BIN)
	rm -rf output
	rm -f energy_four_panel.png
	rm -f energy_final_logplot.png

# -----------------------------------------
# Full clean — also remove MD output files
# -----------------------------------------
distclean: clean
	rm -f output/full_positions.pdb
	rm -f output/cell_positions.pdb
	rm -f output/nbl_positions.pdb
	rm -f output/nbl_pthread_positions.pdb
	rm -f output/final_positions.pdb

	rm -f output/full_energies.csv
	rm -f output/cell_energies.csv
	rm -f output/nbl_energies.csv
	rm -f output/nbl_pthread_energies.csv
	rm -f output/final_energies.csv

	rm -f energy_final_logplot.png
	rm -f energy_four_panel.png

	rm -rf output/*

.PHONY: all run clean distclean output compare plot