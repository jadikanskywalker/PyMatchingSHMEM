#!/bin/bash
cd run
rm *.01 *.b8 *.dem
rm out_frames/*

IN="/mnt/DISCL/home/jadhicks/more-bacon-less-threshold/out/circuits/r=3,d=5,p=0.01,noise=uniform,c=bacon_shor_xx_surgery,q=50,b=X,g=all.stim"

stim analyze_errors \
    --decompose_errors \
    --fold_loops \
    --in $IN \
    > error_model.dem

stim detect \
    --in $IN \
    --shots 100 \
    --obs_out actual_obs_flips.01 \
    --obs_out_format 01 \
    --out detection_events.b8 \
    --out_format b8

echo "Starting serial run..."
start_serial=$(date +%s)
# ~/PyMatchingSHMEM/build/pymatching predict \
#     --dem error_model.dem \
#     --in detection_events.b8 \
#     --in_format b8 \
#     --out predicted_obs_flips.01 \
#     --out_format 01 \
#     > log_serial.out

~/PyMatchingSHMEM/build/pymatching animate \
    --dem_in error_model.dem \
    --dets_in detection_events.b8 \
    --dets_in_format b8 \
    --held_frames_per_event 1 \
    --held_frames_at_start 10 \
    --held_frames_at_end 10 \
    --max_growth_between_frames 25 \
    --pixels_per_unit_length 20 \
    --out_dir out_frames
end_serial=$(date +%s)
serial_time=$((end_serial - start_serial))
echo "Serial run completed in $serial_time seconds."

echo Serial
echo correct predictions:
paste -d " " predicted_obs_flips.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l
echo