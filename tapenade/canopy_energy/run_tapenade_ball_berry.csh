#!/bin/tcsh -f
#
# run_tapenade_ball_berry.csh
#
# One-command differentiation, build and test of one PET method 6 (Ball-Berry canopy energy balance) step
# with the VEGE_FLUX loop1 replica, with Tapenade 3.16: the head ball_berry_head()
# (tapenade/canopy_energy/ball_berry_head.c) runs ball_berry_method6_physics_step()
# (src/ball_berry_canopy_physics.c, src/noahmp_stomata.c, src/ball_berry_struct_copies.c) with every active
# quantity in one array of independents and one of dependents (ball_berry_head.h), the canopy state at the
# start of the step among the independents.
#
# The tests (ball_berry_ad_test.c) replay the inputs dcfe-noahmp-stack gave Ball-Berry for cat-18560 of the
# 01123000 basin (data/, the whole realization period, written by build/stack_driver --core-path
# --dump-slot-inputs), each hour through ball_berry_method6_core_step(), with the divide's configuration:
#   p  the head against the production step at every hour, bitwise
#   t  the tangent against central differences at the test states (one per regime and every 5000th
#      hour), with a step-size scan, over the dependents each step resolves; kinks, directions at a bound
#      of the domain and directions below finite-difference resolution reported apart
#   a  the reverse code against the tangent: the full Jacobian and the dot-product test
#
# Usage, from anywhere:
#     tapenade/canopy_energy/run_tapenade_ball_berry.csh          (experiments pta)
#     tapenade/canopy_energy/run_tapenade_ball_berry.csh pt       (any letters)
#     tapenade/canopy_energy/run_tapenade_ball_berry.csh -h       (this help)
# Without Tapenade the test is built from the primal code alone: p runs, t reports the central
# differences' self-consistency, a is skipped.
#
# Steps:
#   1. stage the differentiated sources and the headers they include in WORK_DIRECTORY/canopy
#   2. Tapenade: tangent, then reverse (with fewer optimizations if Tapenade fails), with
#      tapenade/BBGeneralLib declaring hypot (retried without it if Tapenade fails)
#   3. repair known Tapenade 3.16 output problems; make the generated headers coexist with the
#      originals (as sft-C, CIESM, NoahMP_snow_ablation and CFE3.1)
#   4. print the generated prototypes; read the derivative slots of the passive arguments and their types
#   5. compile; then link everything with the production code and the test
#   6. run; results in RESULTS_DIRECTORY
#
# Everything generated lives in WORK_DIRECTORY or RESULTS_DIRECTORY (both ignored by git).
#
# FLO 10/28: first version (from CFE3.1's tapenade/cfe_hour/run_tapenade_cfe_hour.csh).  ASCII only.

# ===================== directories =====================
set script_directory = $0:h
if ( "$script_directory" == "$0" ) set script_directory = .
set EXPERIMENT_DIRECTORY = `cd $script_directory && pwd`
set REPOSITORY = `cd $EXPERIMENT_DIRECTORY/../.. && pwd`
set DATA_DIRECTORY = $EXPERIMENT_DIRECTORY/data
set WORK_DIRECTORY = $REPOSITORY/tapenade_work/canopy_energy
set RESULTS_DIRECTORY = $REPOSITORY/output/tapenade_canopy_energy
set TAPENADE_HOME_DEFAULT = /user2/ogden/opt/tapenade/tapenade_3.16

# ===================== settings =====================
set EXPERIMENTS = pta
set CONFIGURATION = $DATA_DIRECTORY/cat-18560_ball_berry.cfg
set CANOPY_ENERGY_DUMP = $DATA_DIRECTORY/stack_01123000_cat-18560_canopy_energy.csv.gz
# the head: name, Tapenade head specification, passive arguments that may get derivative slots
set HEAD_NAMES = ( canopy )
set HEAD_SPECIFICATION_canopy = "ball_berry_head(dependent)/(independent)"
set PASSIVE_ARGUMENTS_canopy = ( passive_parameters passive_step )
set PRODUCTION_CFLAGS = "-std=c99 -O2 -Wall -Wextra -g -ffp-contract=off"
set GENERATED_CFLAGS = "-std=c11 -O2 -Wall -Wextra -g -ffp-contract=off"

# ===================== end of settings =====================
if ( $#argv > 0 ) then
    if ( "$1" == "-h" || "$1" == "--help" ) then
        sed -n '2,/^# ====/p' $0 | sed '$d' | sed 's/^# \{0,1\}//'
        exit 0
    endif
    set EXPERIMENTS = "$1"
endif

set nonomatch
cd $REPOSITORY
set object_directory = $WORK_DIRECTORY/objects
set test_program = $WORK_DIRECTORY/ball_berry_ad_test
set staged_names = ( ball_berry_head.c ball_berry_canopy_physics.c ball_berry_struct_copies.c noahmp_stomata.c )
set staged_sources = ( $EXPERIMENT_DIRECTORY/ball_berry_head.c src/ball_berry_canopy_physics.c \
                       src/ball_berry_struct_copies.c src/noahmp_stomata.c )
# the headers the differentiated code includes (their macros are force-included into the generated code)
set staged_headers = ( include/noahmp_stomata.h include/ball_berry_helpers.h include/ball_berry_state.h \
                       include/ball_berry_canopy_physics.h include/ball_berry_struct_copies.h \
                       $EXPERIMENT_DIRECTORY/ball_berry_head.h )
# the production sources: every source of the library (src/*.c without the drivers), and the head
set production_sources = ( )
foreach source_file ( src/*.c )
    if ( "$source_file:t" =~ main_* ) continue
    set production_sources = ( $production_sources $source_file )
end
set production_sources = ( $production_sources $EXPERIMENT_DIRECTORY/ball_berry_head.c )

if ( $?TAPENADE_HOME ) then
    set tapenade_home = $TAPENADE_HOME
else
    set tapenade_home = $TAPENADE_HOME_DEFAULT
endif
set tapenade = $tapenade_home/bin/tapenade
set kit_directory = $tapenade_home/ADFirstAidKit
set have_tapenade = 1
if ( ! -x $tapenade ) then
    echo "NOTE: Tapenade not found at $tapenade (set TAPENADE_HOME): building the primal-only test"
    set have_tapenade = 0
else if ( ! -f $kit_directory/adStack.c || ! -f $kit_directory/adStack.h ) then
    echo "ERROR: reverse mode needs $kit_directory/adStack.c and adStack.h (ADFirstAidKit)"
    exit 1
endif

rm -rf $WORK_DIRECTORY
mkdir -p $object_directory $RESULTS_DIRECTORY
gzip -d -c $CANOPY_ENERGY_DUMP >! $WORK_DIRECTORY/canopy_energy.csv
if ( $status != 0 ) then
    echo "ERROR: cannot read $CANOPY_ENERGY_DUMP"
    exit 1
endif

set objects = ""
set generated_tarball_list = ""
if ( $have_tapenade ) then
    foreach head ( $HEAD_NAMES )
        eval set head_specification = '"$HEAD_SPECIFICATION_'$head'"'
        eval set passive_arguments = '( $PASSIVE_ARGUMENTS_'$head' )'
        set head_function = `echo "$head_specification" | sed 's/(.*//'`
        set head_directory = $WORK_DIRECTORY/$head
        set generated_d = $head_directory/generated_d
        set generated_b = $head_directory/generated_b
        mkdir -p $generated_d $generated_b
        echo "==== ${head}: 1. staging the differentiated sources in $head_directory"
        cp $staged_headers $head_directory/
        foreach source_file ( $staged_sources )
            cp $source_file $head_directory/$source_file:t
        end

        echo "==== ${head}: 2. running Tapenade (tangent mode, then reverse mode)"
        set library = $REPOSITORY/tapenade/BBGeneralLib
        cd $head_directory
        $tapenade -tangent -ext $library -head "$head_specification" -O generated_d $staged_names >& generated_d/tapenade.log
        if ( ! -f generated_d/ball_berry_head_d.c ) then
            echo "NOTE: tapenade -tangent with BBGeneralLib failed (log saved); retrying without it"
            cp generated_d/tapenade.log failed_tangent_with_library.log
            set library = ""
            rm -f generated_d/*
            $tapenade -tangent -head "$head_specification" -O generated_d $staged_names >& generated_d/tapenade.log
        endif
        set reverse_options_used = "none"
        foreach reverse_options ( "" "-nooptim recomputeintermediates" "-nooptim recomputeintermediates -nooptim tbr" )
            rm -f generated_b/*
            if ( "$library" != "" ) then
                $tapenade -reverse $reverse_options -ext $library -head "$head_specification" -O generated_b \
                    $staged_names >& generated_b/tapenade.log
            else
                $tapenade -reverse $reverse_options -head "$head_specification" -O generated_b $staged_names \
                    >& generated_b/tapenade.log
            endif
            if ( -f generated_b/ball_berry_head_b.c ) then
                set reverse_options_used = "$reverse_options"
                if ( "$reverse_options_used" == "" ) set reverse_options_used = "(defaults)"
                break
            endif
            echo "NOTE: tapenade -reverse $reverse_options failed; log saved"
            cp generated_b/tapenade.log failed_reverse_`echo "$reverse_options" | tr -c 'a-z' '_'`.log
        end
        cd $REPOSITORY
        echo "reverse-mode options that worked: $reverse_options_used"
        foreach mode ( d b )
            set generated = $head_directory/generated_$mode
            if ( ! -f $generated/ball_berry_head_$mode.c ) then
                echo "ERROR: Tapenade did not produce $generated/ball_berry_head_$mode.c.  Log:"
                cat $generated/tapenade.log
                exit 1
            endif
            echo "---- Tapenade messages, ${head}, mode $mode (JAVA_TOOL_OPTIONS and file-creation lines omitted):"
            grep -v "JAVA_TOOL_OPTIONS" $generated/tapenade.log | grep -v "^@@ Creat" | grep -v "empty type for result"
        end
        set generated_tarball_list = "$generated_tarball_list $head/generated_d $head/generated_b"

        echo "==== ${head}: 3. repairing generated code, and making the generated headers coexist with the originals"
        foreach generated_source ( $generated_d/*.c $generated_b/*.c )
            set empty_operand_count = `grep -c '?  :' $generated_source`
            if ( $empty_operand_count != 0 ) then
                sed -i 's/?  :/? 0.0 :/g' $generated_source
                echo "   ${generated_source}: supplied 0.0 for $empty_operand_count empty conditional operands"
            endif
        end
        # Tapenade 3.16 declares the dummy derivatives of constant arguments as float where the derivative
        # functions take double * (CIESM's minimum_of_two and maximum_of_two, nordin 10/28). Nothing reads
        # them: declare them double
        foreach generated_source ( $generated_d/*.c $generated_b/*.c )
            set float_dummy_count = `grep -c -E 'float[ ]+dummydiff' $generated_source`
            if ( $float_dummy_count != 0 ) then
                sed -i -E 's/float([ ]+dummydiff)/double\1/g' $generated_source
                echo "   ${generated_source}: $float_dummy_count float dummy derivatives declared double"
            endif
        end
        foreach generated ( $generated_d $generated_b )
            foreach generated_header ( $generated/*_[db].h )
                set base = `basename $generated_header .h | sed 's/_[db]$//'`
                if ( -f $head_directory/$base.h ) then
                    set guard = `grep -m1 '^#ifndef' $head_directory/$base.h | awk '{print $2}'`
                    sed -i "0,/^#define [A-Za-z_]*_[DB]_LOADED"'$'"/s//&\n#define $guard/" $generated_header
                endif
            end
            # (one-line macros only)
            ( echo "/* macros of the original headers, for compiling Tapenade-generated code */" ; \
              cat $staged_headers | grep '^[ 	]*#define [A-Za-z_0-9]*[ 	(][ 	]*[^ 	]' | \
              grep -v '_H$' | grep -v '\\$' | sed 's/^[ 	]*//' | sed 's|[ 	]*/[*].*$||' | awk '\!seen[$2]++' ) >! \
                $generated/original_header_macros.h
        end

        echo "==== ${head}: 4. generated prototypes, and the derivative slots of the passive arguments"
        set head_upper = `echo $head | tr a-z A-Z`
        foreach mode ( d b )
            set generated = $head_directory/generated_$mode
            set mode_upper = `echo $mode | tr a-z A-Z`
            set prototype_file = $head_directory/prototype_$mode.txt
            grep -A24 "void ${head_function}_$mode *(" $generated/ball_berry_head_$mode.c | sed -n '1,/)/p' | \
                tr '\n' ' ' >! $prototype_file
            cat $prototype_file
            echo ""
            set slot_flags = ""
            foreach argument ( $passive_arguments )
                set argument_upper = `echo $argument | tr a-z A-Z`
                # (the pattern is built in pieces: tcsh would read $mode[...] as a word selector)
                set slot_name = "${argument}${mode}"
                set pointer_count = `grep -c -E '[*][ ]*'"$slot_name"'[ ,)]' $prototype_file`
                set any_count = `grep -c -E "$slot_name"'[ ,)]' $prototype_file`
                if ( "$pointer_count" != "0" ) then
                    # the slot's type: Tapenade 3.16 names a structure's derivative type <name>_diff, but gives
                    # the slot the original type when the structure has no differentiable field of its own
                    # (CFE3.1: cfe_options_struct, cfe_volbal_struct)
                    set slot_type = `grep -o -E '[A-Za-z_][A-Za-z_0-9]*[ ]*[*][ ]*'"$slot_name"'[ ,)]' $prototype_file | head -1 | sed 's/[ *].*//'`
                    set slot_flags = "$slot_flags -DBB_HEAD_${mode_upper}_${argument_upper}_SLOT_BY_POINTER"
                    set slot_flags = "$slot_flags -DBB_HEAD_${mode_upper}_${argument_upper}_SLOT_TYPE=$slot_type"
                else if ( "$any_count" != "0" ) then
                    set slot_flags = "$slot_flags -DBB_HEAD_${mode_upper}_${argument_upper}_SLOT_BY_VALUE"
                endif
            end
            if ( "$mode" == "d" ) set adapter = $EXPERIMENT_DIRECTORY/ball_berry_tangent_adapter.c
            if ( "$mode" == "b" ) set adapter = $EXPERIMENT_DIRECTORY/ball_berry_adjoint_adapter.c

            echo "==== ${head}: 5. compiling, mode $mode (slot flags:$slot_flags)"
            foreach generated_source ( $generated/*_$mode.c )
                if ( "$generated_source:t" == "ball_berry_head_$mode.c" ) continue
                set object = $object_directory/${head}_$generated_source:t:r.o
                cc $GENERATED_CFLAGS -include $generated/original_header_macros.h -I$generated -I$head_directory \
                    -I$kit_directory -c $generated_source -o $object
                if ( $status != 0 ) then
                    echo "ERROR: compiling $generated_source failed"
                    exit 1
                endif
                set objects = "$objects $object"
            end
            cc $GENERATED_CFLAGS $slot_flags -include $generated/original_header_macros.h -I$generated -I$head_directory \
                -I$EXPERIMENT_DIRECTORY -I$kit_directory -c $adapter -o $object_directory/$adapter:t:r.o
            if ( $status != 0 ) then
                echo "ERROR: compiling $adapter:t failed (the generated signature differs from the expected one; see the"
                echo "       prototype printed above)"
                exit 1
            endif
            set objects = "$objects $object_directory/$adapter:t:r.o"
        end
    end
    tar czf $RESULTS_DIRECTORY/generated_code.tar.gz -C $WORK_DIRECTORY $generated_tarball_list
    echo "generated code (repaired): $RESULTS_DIRECTORY/generated_code.tar.gz"
    cc -std=c11 -O2 -D_POSIX_C_SOURCE=200809L -I$kit_directory -c $kit_directory/adStack.c -o $object_directory/adStack.o
    if ( $status != 0 ) then
        echo "ERROR: compiling adStack.c failed"
        exit 1
    endif
    set objects = "$objects $object_directory/adStack.o"
endif

echo ""
echo "==== compiling the production code and the test"
set primal_flag = ""
if ( ! $have_tapenade ) set primal_flag = "-DBB_AD_PRIMAL_ONLY"
set production_objects = ""
foreach source_file ( $production_sources $EXPERIMENT_DIRECTORY/ball_berry_ad_test.c )
    set object = $object_directory/$source_file:t:r.o
    cc $PRODUCTION_CFLAGS $primal_flag -Iinclude -I$EXPERIMENT_DIRECTORY -c $source_file -o $object
    if ( $status != 0 ) then
        echo "ERROR: compiling $source_file failed"
        exit 1
    endif
    set production_objects = "$production_objects $object"
end
# the production objects first: the generated files carry copies of some production functions, and
# with duplicate definitions allowed the linker keeps the first
cc -o $test_program $production_objects $objects -lm -Wl,--allow-multiple-definition
if ( $status != 0 ) then
    echo "ERROR: link failed"
    exit 1
endif

echo ""
echo "==== 6. running experiments $EXPERIMENTS"
$test_program $CONFIGURATION $WORK_DIRECTORY/canopy_energy.csv $RESULTS_DIRECTORY $EXPERIMENTS |& \
    tee $RESULTS_DIRECTORY/results_$EXPERIMENTS.txt
echo ""
echo "Results: $RESULTS_DIRECTORY/results_$EXPERIMENTS.txt"
