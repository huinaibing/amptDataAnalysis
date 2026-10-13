# AMPT Data Analysis

ROOT macros for producing O2-compatible flow and mean-pT outputs from AMPT
particle trees. Run the macros from `script/` with ROOT 6 and C++17 support.
The required O2Physics sources are included in `external/o2physics/`; an
O2Physics installation is no longer required. ROOT interprets these sources
directly, so run macros without the ACLiC `+` suffix.

The JSON parser in `dataFrame/json.hpp` is the unmodified nlohmann/json
v3.12.0 `single_include/nlohmann/json.hpp` from the upstream release tag.
It includes all JSON implementation headers, so no system nlohmann/json
installation is required. Its MIT license is in `dataFrame/JSON_LICENSE.MIT`.

## Entry points

- `script/calculate_v2ptrho.cpp` produces the `processData`-style v2-pT-rho
  FlowContainer and `meanptCentNbs` output.
- `script/calculate_c22deltapt.cpp` produces the
  `processDataC22DeltaPt`-style profiles under `c22DeltaPt/`.
- `script/calculate_pidptcorrelations.cpp` produces the
  `processPidPtCorrelations`-style PID mean-pT `FlowContainer`.
- `script/calculate_pid_charged_sc.cpp` produces the pure PID and charged
  two-track moments and the mixed two-PID/two-charged four-track moment for
  SC(v2,PID squared, v2,charged squared), under `pid-charged-sc/`.

These functions accept:

1. input-file JSON path;
2. output ROOT path;
3. maximum files per JSON entry (`-1` means all);
4. maximum JSON entries (`-1` means all);
5. analysis JSON path.

The c22-delta-pT macro has a `usePure` argument before the analysis JSON path.
Its default `-1` reads `use_pure` from JSON; `0` selects PID POI-ref and `1`
selects PID POI-POI.

Quick one-file examples:

```bash
root -l -b -q 'calculate_v2ptrho.cpp("../config/cent_cfg.json", "/tmp/v2ptrho.root", 1, 1, "../config/calculate_v2ptrho.json")'
root -l -b -q 'calculate_c22deltapt.cpp("../config/cent_cfg.json", "/tmp/c22.root", 1, 1, -1, "../config/calculate_c22deltapt.json")'
root -l -b -q 'calculate_pidptcorrelations.cpp("../config/cent_cfg.json", "/tmp/pidpt.root", 1, 1, "../config/calculate_pidptcorrelations.json")'
```

To run the analysis across independent ROOT processes, load ROOT and run
from `script/`:

```bash
./run_parallel.sh calculate_v2ptrho 4 ../config/cent_cfg.json /tmp/v2ptrho.root -1 -1 ../config/calculate_v2ptrho.json
```

The second argument is the number of processes. The same command accepts
`calculate_c22deltapt`, `calculate_pidptcorrelations`, or
`calculate_flowGfwOmegaXi`, `calculate_flowPbpbPikp`, or
`calculate_ptspectra`, `calculate_pid_v2_rp_gap_variance`, or
`calculate_pid_charged_sc`. For `calculate_c22deltapt`, the mode is read from
`c22_delta_pt_output.use_pure` in the analysis JSON: `false` selects PID
POI-ref and `true` selects PID POI-POI Pure. The shell does not accept a
separate mode argument. Each process reads a disjoint
set of input file numbers; `merge_results.cpp` loads the bundled
`FlowContainer` and merges their output after all succeed. The
script prints the location of per-process ROOT files and logs. Direct macro
calls still work; their optional final arguments are `shardIndex` and
`shardCount`.

On the CERN server, generate the input JSON with `update_ampt_json_eos.sh`
from the repository root, after setting its path variables. Then run, for
example:

```bash
cd script
root -l -b -q 'calculate_v2ptrho.cpp("../config/cent_cfg_ampt_eos.json", "/tmp/v2ptrho.root", 1, 1, "../config/calculate_v2ptrho.json")'
./run_parallel.sh calculate_v2ptrho 4 ../config/cent_cfg_ampt_eos.json /tmp/v2ptrho_all.root -1 -1 ../config/calculate_v2ptrho.json
```

The JSON `path` may be an EOS mount path or an XRootD prefix such as
`root://eosuser.cern.ch//eos/user/q/qiuyu/ampt/AMPT_NEW/Result/ampt_<id>_`.
The latter requires ROOT's XRootD support and EOS authentication. JSON
configuration files and output paths in these examples remain local.

For hydro files named `hydro_<id>_<number>.root`, run
`./update_hydro_json_eos.sh` from the repository root on lxplus. Its path
variables default to `/eos/user/q/qiuyu/hydro/result/PbPb5360/particles/`,
and it writes `config/cent_cfg_hydro_eos.json` without changing the local
hydro input list. Like the AMPT version, it uses `eos ls -F` instead of
scanning the EOS mount with `find`.

## CERN HTCondor

For large server runs, use HTCondor so the ROOT processes run on batch
workers. Edit the settings in `config/condor/submit.conf`,
then submit without command-line arguments:

```bash
cd script
./submit_condor.sh
```

The default configuration is located relative to the repository, so the
script can also be run from the repository root as
`./script/submit_condor.sh`. The configuration uses Bash variable syntax;
`REPO_DIR` and `SCRIPT_DIR` are supplied by the script. Use these variables
or absolute paths for file and directory settings.

To keep separate submission settings for another dataset or task, copy
`config/condor/submit.conf`, edit the copy, and pass its path:

```bash
./submit_condor.sh ../config/condor/submit_hydro.conf
```

Set `INPUT_JSON` to the EOS input list generated by
`update_ampt_json_eos.sh` and `EOS_OUTPUT_DIR` to your own writable EOS XRootD
directory. `ANALYSIS_JSON` defaults to the file for `MACRO`; set it explicitly
only to use a custom configuration. The script
converts `/eos/` input prefixes to `root://eosuser.cern.ch//eos/...` in the
packaged copy. It preserves every JSON entry and its original order, so
file sharding and deterministic bootstrap identities match `run_parallel.sh`.
The repository's original JSON and analysis macros are unchanged.
For `calculate_c22deltapt`, set `c22_delta_pt_output.use_pure` in
`ANALYSIS_JSON`; the Condor shells do not define or override this mode.

To submit hydro data with the same analysis, set these variables:

```bash
INPUT_JSON="$REPO_DIR/config/cent_cfg_hydro_eos.json"
OUTPUT_NAME="myAnalysisResultV2PtRho_hydro.root"
EOS_OUTPUT_DIR="root://eosuser.cern.ch//eos/user/q/qiuyu/hydro/condor"
```

`MACRO` still selects the calculation (for example `calculate_v2ptrho`),
and `ANALYSIS_JSON` still selects the physics settings. Use a separate
analysis JSON if hydro needs different selections, centrality settings,
or output axes. ROOT setup and the EOS endpoint use the same defaults;
choose shard and job limits for the dataset size. The hydro ROOT files
must expose the same `particles` tree and branches expected by `dataLoader.h`.

`ROOT_SETUP` selects a fixed ordinary ROOT environment available to the
workers; the default is the CERN CVMFS LCG_108 EL9/gcc13 view (ROOT 6.36.02).
`run.sub` selects EL9/x86_64 workers and uses `JobFlavour = "longlunch"`
(2 hours), with no explicit CPU or memory request. Both analysis and
merge use this submit description. Increase `SHARD_COUNT` if an analysis
shard exceeds the queue's time limit. The merge must also finish
within that limit; otherwise change the flavour before resubmitting.

The submit script packages code, bundled dependencies, and both JSONs.
Each `condor_worker.sh` job unpacks this snapshot into worker scratch and
runs one ROOT process with an explicit shard index. All shards must
succeed, including successful EOS output upload, before DAGMan starts
the separate `merge_results.cpp` job. Each run has its own EOS directory
containing `part_0.root`, ..., and `OUTPUT_NAME`. Failed jobs do not upload
partial analysis results. Nonexistent raw input files are skipped; the
input-count check includes both loaded and skipped files. Other ROOT
errors, input-count mismatches, invalid output files, and EOS upload
failures cause a nonzero exit. DAGMan retries failed nodes up to `RETRIES`
times. After validation, the worker creates the EOS run directory with
`xrdfs mkdir -p` and uploads the ROOT file with `xrdcp --force`, so DAGMan
only considers the node successful after the upload completes.

Keep `SUBMIT_DIR` on AFS (the default is `script/condor_runs/` when the
repository is in AFS), and have a valid Kerberos ticket before submission.
The CERN batch service supplies worker credentials; the submit file does
not inherit the login shell environment. Submit files, DAG state, stdout,
stderr, and scheduling logs remain under `SUBMIT_DIR`. `ON_EXIT` returns
stdout/stderr even when the worker fails, while `transfer_output_files =
""` prevents Condor from requiring a ROOT file for a failed job. CERN does
not permit `stream_output` or `stream_error`. Original data and ROOT results
reside on EOS. `MAX_JOBS` limits queued/running DAG nodes,
independently of the total `SHARD_COUNT`.

Existing run directories contain their own copied submit file, worker,
and code archive; repository changes do not update them. To recover a
failed DAG, copy the updated `run.sub` and `condor_worker.sh` into that run
directory, then run `condor_submit_dag -dorescuefrom 1 analysis.dag` for
`analysis.dag.rescue001`. This keeps completed nodes, reruns failed shards,
then submits `MERGE`. Do not use `-force` or rerun `submit_condor.sh` for
this recovery, as that would start the completed shards again.

For a remote pilot, first set `SHARD_COUNT=1`, `MAX_FILES=1`, and
`MAX_CONFIGS=1`. Verify EOS reads, output transfer, and the automatic merge
in that run before restoring the full input limits. To inspect the
generated archive, submit file and DAG before submission, set
`SUBMIT=false`; the script prints the manual submission command. Monitor
with `condor_q`, inspect the run's `logs/`, and remove the entire DAG with
`condor_rm <DAGMan-job-id>` if needed.

For postprocessing with plain ROOT, load
`external/o2physics/FlowContainerRoot.h` before reading `FlowContainer`
objects. The persisted class names, versions, and output paths match O2.
See `external/o2physics/README.md` for upstream revisions and local changes.

## PID–charged SC inputs

Run from `script/` after loading ROOT (no ACLiC compilation):

```bash
root -n -l -b -q 'calculate_pid_charged_sc.cpp("../config/cent_cfg_ampt.json", "/tmp/pid_charged_sc_moments.root", 1, 1, "../config/calculate_pid_charged_sc.json")'
./run_parallel.sh calculate_pid_charged_sc 4 ../config/cent_cfg_ampt.json /tmp/pid_charged_sc_all.root -1 -1 ../config/calculate_pid_charged_sc.json
```

For Condor, set `MACRO="calculate_pid_charged_sc"` and
`OUTPUT_NAME="pid_charged_sc_moments.root"` in `config/condor/submit.conf`, and choose
the input list and EOS output directory as usual. Its task directory mapping
is `pid-charged-sc`; the generic worker already supports this macro's argument
order and completion message.

The `pid_charged_sc_output` section in
`config/calculate_pid_charged_sc.json` controls eta
subevents, strict cut boundaries, centrality and subsample axes, and event
weights. The independent pT windows reuse `selection.pt_ranges` for charged,
pion, kaon and proton. Defaults select `-0.8 < eta < -0.4` (A) and
`0.4 < eta < 0.8` (B), 1% centrality bins and 30 disjoint subsamples.

Each `pid-charged-sc/{Pion,Kaon,Proton}/` directory contains:

- `pTwoPID`: PID A–B two-track moment;
- `pTwoCharged`: charged A–B moment on this species' valid event sample;
- `pFourPIDCharged`: distinct PID-A, charged-A, PID-B, charged-B moment
  with harmonics `+2,+2,-2,-2`;
- the corresponding `*ByGroup` TProfile2D objects (centrality, subsample);
- `hValidEvents` and `hSumPIDPairs/hSumChargedPairs/hSumMixedQuadruplets`.

Only events with at least one valid mixed quadruplet contribute to all three
moments. Correlations are evaluated with the bundled GFW cumulant framework.
`PidChargedScCorrConfigManager` in `corrConfigManager.h` configures the charged,
PID and overlap regions. For example, the pion mixed four-track config is
`scPoiPiA scRefA | scOlPiA {2 2} scPoiPiB scRefB | scOlPiB {-2 -2}`.
GFW recursively subtracts repeated tracks using the actual selection overlap
on each side. `Calculate(config, 0, false)` supplies the numerator and
`Calculate(config, 0, true)` the distinct tuple count. All four tracks are
distinct, including when PID is included in the charged sample. No
event-by-event square root or POI-reference ratio is used.

With the default `event_weight="combinations"`, each moment uses its own
distinct pair/quadruplet count, as the other flow macros do. The optional
`"mixed_quadruplets"` mode weights all three moments by the mixed quadruplet
count, for a covariance with one common multiplicity-weighted event measure.
Neither mode uses unit event weights. Do not mix modes when merging results.

Downstream, form SC = `pFourPIDCharged - pTwoPID*pTwoCharged` in the original
narrow centrality bins, then combine the results with an explicit weight
convention. NSC divides SC by `pTwoPID*pTwoCharged` when that denominator is
valid. Use the same operation for each subsample to obtain statistical
uncertainties. With native tuple weights, the three averages have different
event weights: this is the conventional correlator-weighted SC estimator,
rather than an exact covariance under one common event probability measure.
The two-subevent moments also retain any longitudinal flow decorrelation.

`hEventCount/hCent` are unweighted event-count QA histograms. Fixed `TNamed`
metadata record the GFW engine, weight mode, harmonic, layout, event sample, subsample
count and analysis JSON. Only mergeable upstream components are saved; SC,
NSC and their uncertainties are left to downstream processing.

The numeric regression macro `test/test_pid_charged_sc.cpp` independently
enumerates all distinct tuples, checks partial overlaps and unequal pT cuts,
and verifies both weight modes and subsample sums. Run it from `test/` with
`root -n -l -b -q 'test_pid_charged_sc.cpp("../config/calculate_pid_charged_sc.json")'`.

## Configuration

Configuration files live in `config/`. Input-file lists such as
`cent_cfg_ampt.json` and `cent_cfg_hydro.json` remain separate from analysis
settings. `config/condor/submit.conf` holds Condor submission settings.
Each analysis macro defaults to its own file:

| Macro | Analysis configuration | Output section |
| --- | --- | --- |
| `calculate_v2ptrho` | `config/calculate_v2ptrho.json` | `v2_pt_rho_output` |
| `calculate_c22deltapt` | `config/calculate_c22deltapt.json` | `c22_delta_pt_output` |
| `calculate_pidptcorrelations` | `config/calculate_pidptcorrelations.json` | `pid_pt_correlations_output` |
| `calculate_flowGfwOmegaXi` | `config/calculate_flowGfwOmegaXi.json` | Fixed O2 axes in the macro |
| `calculate_flowPbpbPikp` | `config/calculate_flowPbpbPikp.json` | `flowPbpbPikp_output` |
| `calculate_ptspectra` | `config/calculate_ptspectra.json` | `pt_spectra_output` |
| `calculate_pid_v2_rp_gap_variance` | `config/calculate_pid_v2_rp_gap_variance.json` | `pid_v2_rp_gap_variance_output` |
| `calculate_pid_charged_sc` | `config/calculate_pid_charged_sc.json` | `pid_charged_sc_output` |

Each file independently contains `selection`,
`centrality_from_impact_parameter`, and `bootstrap`, plus its own output
section where applicable. Edit the file for the task you are running;
changes to pT ranges, eta cuts, centrality conversion or bootstrap seed do
not affect other tasks. Existing parameter values were preserved when the
shared configuration was split.

Direct macro calls can still pass an explicit analysis JSON path.
`run_parallel.sh` defaults to `../config/<macro>.json`, and
`submit_condor.sh` defaults to `$REPO_DIR/config/$MACRO.json`. Both retain
support for a custom analysis JSON. Each macro with configurable output
requires its own output section, while unrelated output sections are optional.

`v2_pt_rho_output.qa_axes` controls `hCent/hPt/hPhi/hEta` independently of
the analysis axes. `pid_pt_correlations_output` controls strict pT/eta
boundaries, centrality bins, and bootstrap subsamples for
`calculate_pidptcorrelations.cpp`. Every axis uses one of these two forms:

```json
{"binning": "uniform", "bins": 300, "min": 0.0, "max": 3.0}
{"binning": "variable", "edges": [0.0, 0.2, 0.5, 1.0, 3.0]}
```

## Structure

- `analysisConfig.h`: JSON configuration reader, validation, and PID metadata.
- `analysisUtils.h`: shared event selection, GFW filling, correlation results,
  centrality conversion, and bounded input iteration.
- `corrConfigManager.h`: GFW regions and named correlation configurations.
- `dataLoader.h`: streaming ROOT `particles` row reader.
- `eventManager.h`: groups streamed rows into one in-memory event at a time.
- `dataFrame/`: plain data types and JSON configuration reader.
- `utils.h`: FlowContainer and v2-pT-rho output helpers.
- `pidChargedSc.h`: PID/charged GFW filling, correlator evaluation and mergeable profiles.
- `selection.h`: compatibility wrappers for older macros.
- `external/o2physics/`: bundled ROOT-only flow dependencies and license.

The flow eta gap and mean-pT eta interval are independent. Species-specific
pT ranges are also independent, even where their current defaults are equal.
