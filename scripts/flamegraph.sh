#!/bin/bash

FG_DIR=~/sw/el9-x86_64/FlameGraph

IN=~/PyMatchingSHMEM/logs/vtune_report_topdown.csv
OUT=~/PyMatchingSHMEM/logs/vtune_report_topdown.svg

sed -i '1{/^war:/d}' $IN

perl $FG_DIR/stackcollapse-vtune.pl $IN | perl $FG_DIR/flamegraph.pl --width 1400 --fontsize 10 > $OUT
