import os

# Must be set before the extension's OpenMP runtime spins up threads (inside
# main_predict), which caps num_threads at min(OMP_NUM_THREADS, num_partitions).
os.environ["OMP_NUM_THREADS"] = "1"

import stim
import pymatching as pm

# Equivalent of:
#   stim --rounds 72 --distance 9 --after_clifford_depolarization 0.01 \
#        --code surface_code --task rotated_memory_x > circuit.stim
circuit = stim.Circuit.generated(
    "surface_code:rotated_memory_x",
    rounds=72,
    distance=9,
    after_clifford_depolarization=0.01,
)
circuit.to_file("circuit.stim")

# Equivalent of:
#   stim analyze_errors --decompose_errors --fold_loops --in circuit.stim > error_model.dem
dem = circuit.detector_error_model(decompose_errors=True, flatten_loops=False)
dem.to_file("error_model.dem")

# Equivalent of:
#   stim detect --in circuit.stim --shots 1000 \
#        --obs_out actual_obs_flips.01 --obs_out_format 01 \
#        --out detection_events.b8 --out_format b8
sampler = circuit.compile_detector_sampler()
detection_events, obs_flips = sampler.sample(shots=1000, separate_observables=True)

stim.write_shot_data_file(
    data=detection_events,
    path="detection_events.b8",
    format="b8",
    num_detectors=circuit.num_detectors,
)
stim.write_shot_data_file(
    data=obs_flips,
    path="actual_obs_flips.01",
    format="01",
    num_observables=circuit.num_observables,
)

# Decode with the threaded partitioned decoder: 4 threads (via OMP_NUM_THREADS
# above), 9 rounds per partition.
pm._cpp_pymatching.main(
    command_line_args=[
        "predict",
        "--dem", "error_model.dem",
        "--in", "detection_events.b8",
        "--in_format", "b8",
        "--out", "predicted_obs_flips.01",
        "--out_format", "01"
    ]
)
