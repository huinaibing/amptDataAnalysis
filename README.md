# AMPT Data Analysis

ROOT macros for producing O2-compatible flow and mean-pT outputs from AMPT
particle trees. Run the macros from `script/` inside an O2Physics environment.

## Entry points

- `script/calculate_v2ptrho.cpp` produces the `processData`-style v2-pT-rho
  FlowContainer and `meanptCentNbs` output.
- `script/calculate_c22deltapt.cpp` produces the
  `processDataC22DeltaPt`-style profiles under `c22DeltaPt/`.

Both functions accept:

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
```

## Configuration

Configuration files live in `config/`:

- `config/cent_cfg.json` lists AMPT input prefixes and file counts.
- `config/config.json` contains analysis selections and output axes.

`config.json` controls the flow subevent eta boundaries, the independent
mean-pT eta interval, charged/pion/kaon/proton pT ranges, strict or inclusive
cut boundaries, impact-parameter centrality conversion, bootstrap seed,
POI-POI versus POI-ref mode, and the centrality/mean-pT/bootstrap axes for both
output layouts. Every axis uses one of these two forms:

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

The flow eta gap and mean-pT eta interval are independent. Species-specific
pT ranges are also independent, even where their current defaults are equal.
