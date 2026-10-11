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

All three functions accept:

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
root -l -b -q 'calculate_v2ptrho.cpp("../config/cent_cfg.json", "/tmp/v2ptrho.root", 1, 1, "../config/config.json")'
root -l -b -q 'calculate_c22deltapt.cpp("../config/cent_cfg.json", "/tmp/c22.root", 1, 1, -1, "../config/config.json")'
root -l -b -q 'calculate_pidptcorrelations.cpp("../config/cent_cfg.json", "/tmp/pidpt.root", 1, 1, "../config/config.json")'
```

To run the analysis across independent ROOT processes, load ROOT and run
from `script/`:

```bash
./run_parallel.sh calculate_v2ptrho 4 ../config/cent_cfg.json /tmp/v2ptrho.root -1 -1 ../config/config.json
```

The second argument is the number of processes. The same command accepts
`calculate_c22deltapt`, `calculate_pidptcorrelations`, or
`calculate_flowGfwOmegaXi`, `calculate_flowPbpbPikp`, or
`calculate_ptspectra`. For `calculate_c22deltapt`, an optional final
argument selects `usePure` (`-1`, `0`, or `1`). Each process reads a disjoint
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
root -l -b -q 'calculate_v2ptrho.cpp("../config/cent_cfg_ampt_eos.json", "/tmp/v2ptrho.root", 1, 1, "../config/config.json")'
./run_parallel.sh calculate_v2ptrho 4 ../config/cent_cfg_ampt_eos.json /tmp/v2ptrho_all.root -1 -1 ../config/config.json
```

The JSON `path` may be an EOS mount path or an XRootD prefix such as
`root://eosuser.cern.ch//eos/user/q/qiuyu/ampt/AMPT_NEW/Result/ampt_<id>_`.
The latter requires ROOT's XRootD support and EOS authentication. JSON
configuration files and output paths in these examples remain local.

For postprocessing with plain ROOT, load
`external/o2physics/FlowContainerRoot.h` before reading `FlowContainer`
objects. The persisted class names, versions, and output paths match O2.
See `external/o2physics/README.md` for upstream revisions and local changes.

## Configuration

Configuration files live in `config/`:

- `config/cent_cfg.json` lists AMPT input prefixes and file counts.
- `config/config.json` contains analysis selections and output axes.

`config.json` controls the flow subevent eta boundaries, the independent
mean-pT eta interval, charged/pion/kaon/proton pT ranges, strict or inclusive
cut boundaries, impact-parameter centrality conversion, bootstrap seed,
POI-POI versus POI-ref mode, each macro's output axes, and the
`v2_pt_rho_output.qa_axes` used by `hCent/hPt/hPhi/hEta`. The
`pid_pt_correlations_output` section independently controls strict pT/eta
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
- `selection.h`: compatibility wrappers for older macros.
- `external/o2physics/`: bundled ROOT-only flow dependencies and license.

The flow eta gap and mean-pT eta interval are independent. Species-specific
pT ranges are also independent, even where their current defaults are equal.
