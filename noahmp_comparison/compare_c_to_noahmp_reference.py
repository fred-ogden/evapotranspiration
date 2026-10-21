#!/usr/bin/env python3
# ------------------------------------------------------------------------------------------
# compare_c_to_noahmp_reference.py
#
# Called by compare_C_to_Fortran.sh. Compares one C method 6 solution summary with the
# Fortran NWM Noah-MP reference summaries (r4 = 32-bit REAL as NWM, r8 = 64-bit), and
# optionally the C exact-mode iteration trace with the Fortran r4 iteration trace.
#
# Usage:
#   compare_c_to_noahmp_reference.py [--maximum-relative-difference-vs-r8 X]
#                                    c_summary.csv reference_summary_r4.csv reference_summary_r8.csv
#                                    [c_iteration_trace.csv fortran_iteration_trace_r4.csv]
# With --maximum-relative-difference-vs-r8 the exit status is 1 if any compared variable
# differs from the r8 reference by more than X (relative; absolute where the reference is 0).
# ------------------------------------------------------------------------------------------
import csv
import sys

maximum_relative_difference_vs_r8 = None
if len(sys.argv) > 2 and sys.argv[1] == "--maximum-relative-difference-vs-r8":
    maximum_relative_difference_vs_r8 = float(sys.argv[2])
    del sys.argv[1:3]
if len(sys.argv) not in (4, 6):
    sys.exit("usage: compare_c_to_noahmp_reference.py c_summary.csv reference_r4.csv reference_r8.csv "
             "[c_trace.csv fortran_trace_r4.csv]")

compared_variable_names = ["BTRAN", "FSNO", "SNOWH", "FWET", "LAI", "SAI", "ELAI", "ESAI",
                           "FSUN", "LAISUN", "LAISHA", "PARSUN", "PARSHA", "SAV", "SAG", "FSR",
                           "RSSUN", "RSSHA", "PSNSUN", "PSNSHA", "TV", "TAH", "EAH", "TG", "CM",
                           "IRC", "SHC", "EVC", "TR", "ETRAN_MM_PER_H", "RAHG", "GROUND_LW_ABSORBED"]


def read_name_value_rows(file_path):
    """Return {name: value_text} from 'name,value[,...]' rows; other lines are skipped."""
    name_value_table = {}
    with open(file_path, "r", encoding="ascii") as input_file:
        for line_text in input_file:
            fields = [field.strip() for field in line_text.split(",")]
            if len(fields) >= 2 and fields[0] and fields[0] != "noahmp_name":
                name_value_table[fields[0]] = fields[1]
    return name_value_table


def relative_difference(test_value, reference_value):
    """(test - reference)/|reference|, or the test value itself when the reference is zero."""
    if reference_value == 0.0:
        return test_value
    return (test_value - reference_value) / abs(reference_value)


c_summary = read_name_value_rows(sys.argv[1])
reference_r4 = read_name_value_rows(sys.argv[2])
reference_r8 = read_name_value_rows(sys.argv[3])

print("solver: %s, iterations: %s" % (c_summary.get("SOLVER", "?"), c_summary.get("ITERATIONS", "?")))
print("%-15s %22s %22s %22s %11s %11s" % ("name", "C", "Fortran r4", "Fortran r8", "rel vs r4", "rel vs r8"))
largest_relative_difference_vs_r4 = 0.0
largest_relative_difference_name = ""
largest_relative_difference_vs_r8 = 0.0
largest_relative_difference_vs_r8_name = ""
for variable_name in compared_variable_names:
    if variable_name not in c_summary or variable_name not in reference_r8:
        continue
    c_value = float(c_summary[variable_name])
    r4_value = float(reference_r4[variable_name])
    r8_value = float(reference_r8[variable_name])
    difference_vs_r4 = relative_difference(c_value, r4_value)
    print("%-15s %22.14e %22.14e %22.14e %11.2e %11.2e" % (
        variable_name, c_value, r4_value, r8_value, difference_vs_r4, relative_difference(c_value, r8_value)))
    difference_vs_r8 = relative_difference(c_value, r8_value)
    if abs(difference_vs_r8) > largest_relative_difference_vs_r8:
        largest_relative_difference_vs_r8 = abs(difference_vs_r8)
        largest_relative_difference_vs_r8_name = variable_name
    if abs(difference_vs_r4) > largest_relative_difference_vs_r4:
        largest_relative_difference_vs_r4 = abs(difference_vs_r4)
        largest_relative_difference_name = variable_name
print("largest |relative difference| vs r4: %.2e (%s)" % (largest_relative_difference_vs_r4,
                                                          largest_relative_difference_name))
print("largest |relative difference| vs r8: %.2e (%s)" % (largest_relative_difference_vs_r8,
                                                          largest_relative_difference_vs_r8_name))

if len(sys.argv) == 6:
    with open(sys.argv[4], "r", encoding="ascii") as c_trace_file:
        c_trace_rows = list(csv.DictReader(c_trace_file))
    with open(sys.argv[5], "r", encoding="ascii") as fortran_trace_file:
        fortran_trace_rows = list(csv.DictReader(fortran_trace_file))
    print("")
    print("per-iteration trace, C exact mode vs Fortran r4 (%d and %d iterations)"
          % (len(c_trace_rows), len(fortran_trace_rows)))
    print("%4s %16s %16s %12s %12s %12s %12s" % ("iter", "TV C", "TV Fortran", "RAHC C", "RAHC Fortran",
                                                "TR C", "TR Fortran"))
    for c_row, fortran_row in zip(c_trace_rows, fortran_trace_rows):
        print("%4s %16.9f %16.9f %12.6f %12.6f %12.6f %12.6f" % (
            c_row["iteration"].strip(),
            float(c_row["vegetation_temperature_k"]), float(fortran_row["vegetation_temperature_k"]),
            float(c_row["above_canopy_heat_resistance_s_per_m"]),
            float(fortran_row["above_canopy_heat_resistance_s_per_m"]),
            float(c_row["transpiration_heat_w_per_m2"]), float(fortran_row["transpiration_heat_w_per_m2"])))
    if len(c_trace_rows) != len(fortran_trace_rows):
        print("WARNING: the two codes took different numbers of iterations")

if maximum_relative_difference_vs_r8 is not None:
    if largest_relative_difference_vs_r8 > maximum_relative_difference_vs_r8:
        print("FAIL: largest |relative difference| vs r8 %.2e exceeds %.2e" % (largest_relative_difference_vs_r8,
                                                                            maximum_relative_difference_vs_r8))
        sys.exit(1)
    print("PASS: every compared variable is within %.2e of the r8 reference" % maximum_relative_difference_vs_r8)
