#!/bin/bash

if [[ -n ${1+x} ]]; then
    do_serial=false
else
    do_serial=true
fi

FG_DIR=~/sw/el9-x86_64/FlameGraph

IN=~/PyMatchingSHMEM/run/logs/serial/vtune_report_topdown.csv
OUT=~/PyMatchingSHMEM/run/logs/serial/vtune_report_topdown_serial.svg

INP=~/PyMatchingSHMEM/run/logs/parallel/vtune_report_topdown.csv
OUTP=~/PyMatchingSHMEM/run/logs/parallel/vtune_report_topdown_parallel.svg

# sed -i '1{/^war:/d}' $IN
sed -i '1{/^war:/d}' $INP

if $do_serial; then
    perl $FG_DIR/stackcollapse-vtune.pl $IN | perl $FG_DIR/flamegraph.pl --width 1400 --fontsize 10 > $OUT
fi
perl $FG_DIR/stackcollapse-vtune.pl $INP | perl $FG_DIR/flamegraph.pl --width 1400 --fontsize 10 > $OUTP
