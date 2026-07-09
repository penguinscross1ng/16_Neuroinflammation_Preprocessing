# STAT3-d1mLemon Fiber Photometry Preprocessing

This repository contains MATLAB code for preprocessing and summarizing 24-hour intermittent fiber photometry recordings from mice expressing STAT3-d1mLemon in somatosensory cortex during an LPS-induced inflammation experiment.

The current analysis script is `STAT3_24h.m`.

The file `preprocess_sleep_data.m` is included for review/reference only. It was downloaded from Yue Zhao's `preprocess_sleep_data` repository: <https://github.com/yzhaoinuw/preprocess_sleep_data/blob/main/preprocess_sleep_data.m>. The normalization strategy in `STAT3_24h.m` borrows the reference script's 405-to-465 fitting plus zero-phase smoothing pattern, but this repository does not run the full sleep preprocessing workflow.

## Experiment Outline

- Mice receive local viral injection of STAT3-d1mLemon in somatosensory cortex.
- Fiber photometry recording is performed 3 weeks after surgery.
- At minute 0 of the recording, mice receive intraperitoneal lipopolysaccharide (LPS) to induce an inflammatory response.
- Mice are connected to a TDT fiber photometry system and recorded for 24 hours.
- Recording is intermittent to reduce excessive photobleaching:
  - 96 total cycles
  - 1 cycle = 5 minutes active recording + 10 minutes inactive/off time
  - Total cycle period = 15 minutes
- Current channel mapping in the script:
  - Channel A: `x405A` control and `x465A` signal
  - Channel C: `x405C` control and `x465C` signal

## Analysis Method

The goal is to summarize STAT3-d1mLemon-related fiber photometry signal changes across the 24-hour LPS response while preserving cycle-to-cycle differences.

For each channel and each 5-minute active recording window:

1. Extract the 405 nm control trace and the 465 nm signal trace.
2. Fit the 405 control to the 465 signal using a linear model:

   ```matlab
   fitted405 = slope * f405 + intercept
   ```

3. Smooth the fitted 405 control using a zero-phase moving-average filter.
   - Current smoothing order: `1000` samples
   - Implemented with `filtfilt`
   - This follows the smoothing style used in the referenced `preprocess_sleep_data.m` workflow.

4. Compute percent dF/F:

   ```matlab
   dFF_percent = (f465 - fitted405) ./ fitted405 * 100
   ```

5. Smooth the percent dF/F trace using the same zero-phase moving-average filter.
6. Keep the middle 3 minutes of each 5-minute recording window for summary:
   - Start: 60 seconds after LED-on
   - End: 240 seconds after LED-on

The middle 3-minute window is used to avoid possible edge effects from the beginning and end of each active recording period.

### Z-Score Reference

Because there is no LPS-free baseline recording in this dataset, the first hour is used as an operational within-session reference.

In the script this is:

```matlab
zscore_reference_cycles = 1:4;
```

Those four cycles correspond to the first hour:

- Cycle 1: 0-15 minutes
- Cycle 2: 15-30 minutes
- Cycle 3: 30-45 minutes
- Cycle 4: 45-60 minutes

The z-score reference is computed from the middle 3-minute dF/F samples of cycles 1-4:

```matlab
z_dFF = (dFF_percent - baseline_mean_dFF) / baseline_std_dFF
```

This is a fixed reference for each channel/session. It is not recomputed separately for every 5-minute window, because per-window z-scoring would force each cycle onto its own mean and standard deviation and would reduce interpretability of cycle-to-cycle changes.

Important interpretation note: the first hour occurs after LPS injection, so it should be described as a first-hour reference or operational baseline, not as a true untreated/LPS-free baseline.

## Outputs

The script writes:

- `TDT_dFF_summary.xlsx`
- `TDT_dFF_summary.mat`

The summary table includes, for each channel and cycle:

- Mean/median/max/min/standard deviation of smoothed percent dF/F in the middle 3-minute window
- Mean/median/max/min/standard deviation of z-scored dF/F in the same window
- 405/465/fitted405 middle-window means
- Linear fit slope and intercept
- Smoothing filter order
- Z-score reference mean and standard deviation
- Z-score reference label
- Whether a cycle is part of the z-score reference window

## Limitations And Constraints

- The first hour is used as a practical reference because no LPS-free baseline recording is available. It may already contain early LPS-related biology.
- The script currently fits 405 to 465 separately for each 5-minute active recording window. This may differ from workflows that fit across a full continuous trace or a manually selected fit interval.
- The smoothing order is defined in samples, not seconds. If the sampling rate changes across datasets, the effective smoothing duration changes too.
- Downsampling from the referenced workflow is intentionally not used here.
- The script exports summary tables only; it does not save full normalized traces.
- Artifact rejection is minimal. The script checks for enough samples and finite values, but it does not automatically remove motion artifacts, LED transients, saturation, dropped samples, or abnormal cycles.
- The first and last minute of each 5-minute recording are excluded from summaries, but they are still used in the per-cycle 405-to-465 fit.
- The current script uses local filesystem paths for the TDT SDK, input block, and output directory. These need to be edited for each collaborator's machine.
- Data files are intentionally not committed to this repository.

## Information Still Needed

To make the analysis easier for collaborators to review and reproduce, please add or document:

- Mouse metadata: ID, sex, age, strain, genotype/background, and experimental group.
- Virus details: serotype, promoter, STAT3-d1mLemon construct version, titer, injection volume, injection coordinates, injection depth, and laterality.
- Implant/fiber details: fiber location, ferrule/fiber specs, histology confirmation, and inclusion/exclusion criteria for expression and placement.
- LPS details: dose, concentration, injection volume, lot/vendor, vehicle, and whether injection occurs immediately before the first active recording cycle or during an already-started recording.
- TDT acquisition details: sampling rate, LED powers, excitation wavelengths, demodulation/settings, store names for all animals, and any acquisition filters.
- Experimental context: light/dark cycle, home cage versus other chamber, handling/anesthesia status, food/water access, and behavioral/sleep annotations if available.
- Analysis review decisions: whether 405-to-465 fitting should be per-cycle or fixed across a longer interval, whether smoothing should be defined by seconds instead of samples, and whether full normalized traces should be saved for quality control.
- Statistical plan: animal-level aggregation, treatment/control comparison strategy, multiple-comparison handling, and planned exclusion rules.

## References For Review

The dF/F smoothing approach is based on the local reference script `preprocess_sleep_data.m`, copied from <https://github.com/yzhaoinuw/preprocess_sleep_data/blob/main/preprocess_sleep_data.m>. That script smooths the fitted 405 control and the final percent dF/F trace using a zero-phase moving-average filter. In this repo, the file is committed only so collaborators can inspect the reference method alongside `STAT3_24h.m`.

For biological interpretation of the first-hour reference, collaborators should remember that systemic LPS can induce early brain cytokine responses. Examples include:

- Gabellec et al., 1995: IL-1 family mRNA induction in mouse brain after LPS.
- Nadeau and Rivest, 1999: TNF-alpha gene regulation in brain after systemic immune challenge.
- Laye et al., 2000: hypothalamic cytokine expression measured 1 hour after LPS.
