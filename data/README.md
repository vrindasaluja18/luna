# Datasets

This directory holds public datasets for evaluation. The dataset payloads are
not committed (see `.gitignore`); only this description is tracked.

None of these are wired up yet. The `slamcore` CLI currently reads no external
data, so only the built-in simulator is usable end to end. These directories
are placeholders for the evaluation harness that is still to be written.

| Directory | Dataset | Expected layout |
|---|---|---|
| `euroc/` | EuRoC MAV (ASL) | `mav0/` per sequence, `cam0/`, `imu0/` |
| `kitti/` | KITTI odometry | `sequences/XX/`, `poses/XX.txt`, `calib.txt` |
| `tumvi/` | TUM VI | `dataset-*/mav0/`, `dso/` |
| `malaga/` | MALAGA (parking / urban) | stereo pairs plus GPS/INS ground truth |

To populate, download from the official source and extract so that the layout
above is preserved. Nothing in the build reads these paths yet, so an empty
directory is not an error.
