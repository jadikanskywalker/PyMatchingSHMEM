#!/bin/bash
cd ~/PyMatchingSHMEM/run

rm -r out_frames out_frames_png
rm outpt_video.mp4

stim gen \
    --rounds=8 \
    --distance=8 \
    --after_clifford_depolarization=0.05 \
    --code repetition_code \
    --task memory \
    > circuit.stim
stim analyze_errors \
    --decompose_errors \
    --fold_loops \
    --in circuit.stim \
    > error_model.dem
stim detect \
    --in circuit.stim \
    --shots 1 \
    --obs_out actual_obs_flips.01 \
    --obs_out_format 01 \
    --out detection_events.b8 \
    --out_format b8

../build/pymatching animate \
    --dem_in error_model.dem \
    --dets_in detection_events.b8 \
    --dets_in_format b8 \
    --held_frames_per_event 1 \
    --held_frames_at_start 10 \
    --held_frames_at_end 10 \
    --max_growth_between_frames 25 \
    --pixels_per_unit_length 20 \
    --out_dir out_frames

# mkdir -p out_frames_png
# for file in out_frames/*.svg; do
#     magick "$file" -resize 1024x768 -density 300 "out_frames_png/$(basename "$file" .svg).png" 
# done

# ffmpeg \
#     -framerate 10 \
#     -pattern_type glob \
#     -i 'out_frames_png/*.png' \
#     -vf scale=1024:-1 \
#     -c:v mpeg4 \
#     -q:v 2 \
#     output_video.mp4