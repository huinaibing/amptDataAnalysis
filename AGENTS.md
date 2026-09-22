# AMPT 数据分析项目协作规则

本文件适用于 `/home/huinaibing/git_repo/amptDataAnalysis` 及其子目录。修改文件前先确认文件用途，并尽量保持改动范围最小。上级 `/home/huinaibing/AGENTS.md` 中的通用规则仍然适用。

## 项目结构

- 主函数和入口脚本放在 `script/`。
- 该目录之外的文件通常是头文件或配置文件；除非用户明确要求，不要把主函数放到那里，也不要随意改变这些文件的职责。
- `analysisConfig.h` 负责配置读取和 PID/output metadata，`analysisUtils.h` 负责共享选择、中心度转换、GFW 填充和输入迭代，`corrConfigManager.h` 维护与 O2 GFW 对应的 region 和 correlation，`dataLoader.h`/`eventManager.h` 负责流式读取和按事件分组，`utils.h` 提供 FlowContainer 和 v2-pT-rho 填充工具。公共逻辑应留在这些头文件中，四个主入口只负责各自 workflow 的编排和输出。

## 与 O2 任务的接口契约

`script/` 下的 ROOT 入口宏分别对应 O2 主函数：`pidFlowPtCorr.cxx` 的三个主要数据入口对应 `calculate_v2ptrho.cpp`、`calculate_c22deltapt.cpp`、`calculate_pidptcorrelations.cpp`，`flowGfwOmegaXi.cxx` 的 `processData` 对应 `calculate_flowGfwOmegaXi.cpp`。这些入口使用相同的分析量定义，并把结果写成 O2/AnalysisResults-compatible 的 ROOT 结构。pidFlowPtCorr 对应输出的顶层任务目录统一为 `pid-flow-pt-corr`；flowGfwOmegaXi 对应输出的顶层任务目录为 `flow-gfw-omega-xi`。对应范围内的子目录、FlowContainer 名、profile/histogram 名、轴含义和 bootstrap 组织方式必须与 O2 主函数保持一致。同一个 downstream ROOT 后处理脚本应能按相同路径读取 O2 数据结果和 AMPT 结果。

| AMPT 入口脚本和函数 | 对应的 O2 主函数 | 对应输出 |
| --- | --- | --- |
| `script/calculate_v2ptrho.cpp` / `calculate_v2ptrho(...)` | `processData(...)` | `pid-flow-pt-corr` 下的 `FlowContainerCharged`、`FlowContainerPi`、`FlowContainerKa`、`FlowContainerPr`，以及 `meanptCentNbs` |
| `script/calculate_c22deltapt.cpp` / `calculate_c22deltapt(...)` | `processDataC22DeltaPt(...)` | `pid-flow-pt-corr/hEventCount/processDataC22DeltaPt` 和 `pid-flow-pt-corr/c22DeltaPt` 下的 charged/PID profiles |
| `script/calculate_pidptcorrelations.cpp` / `calculate_pidptcorrelations(...)` | `processPidPtCorrelations(...)` | `pid-flow-pt-corr/FlowContainerPidPtCorr` 和 `pid-flow-pt-corr/hEventCount/processPidPtCorrelations` |
| `script/calculate_flowGfwOmegaXi.cpp` / `calculate_flowGfwOmegaXi(...)` | `processData(...)` | `flow-gfw-omega-xi` 下的 charged flow 基础图（`hPhi/hEta/hPt/hCent/hMult` 等）、`c22/c24/c32/c22Full` 及 charged dpt profiles、Xi/Omega/K0s/Lambda 的 dpt profiles 与不变质量谱 |

修改 O2 主函数的结果对象名、目录、profile 列表、轴顺序、中心度 bin 或 bootstrap 定义时，要检查并同步修改对应的 AMPT 入口；修改 AMPT 输出时也要反向核对 O2 主函数。不要在两边为同一个物理量创建不同名字，也不要让 downstream 脚本通过 AMPT/O2 文件名分支来兼容本可保持一致的结构。

AMPT truth 数据没有 reconstructed `sel8`、detector PID、track-quality、NUA/NUE 和 CCDB 校正。AMPT 入口使用 PDG code 区分带电粒子/粒子种类，以单位粒子权重填充，并通过 impact parameter 转换中心度；`calculate_flowGfwOmegaXi.cpp` 额外用 PDG code 和不变质量窗口重建 K0s、Lambda、Xi、Omega。“两边结构相同”指 ROOT 输出接口和物理量定义相同，不表示 detector-level 选择步骤或 event-count 每个 cut bin 的实际含义完全相同。AMPT 当前只生成适用于 truth 输入的 charged、pion、kaon、proton 及重建 V0/级联粒子结果，不生成 detector PID 失败对应的 unidentified 结果或 O2 QA/校正输出。

## 四个入口脚本

### `calculate_v2ptrho.cpp`

- 对应 `processData(...)`，用于生成 charged 和 pion/kaon/proton 的 flow、mean-pT、pT moments 以及 `v_2-[p_T]` 相关量。
- charged 输出对象名为 `FlowContainerCharged`，PID 输出为 `FlowContainerPi`、`FlowContainerKa`、`FlowContainerPr`。PID 的两个 eta orientation 仍分开填入同一个 profile，以保持 O2 `processData` 的 profile error bookkeeping。
- 同时在 `pid-flow-pt-corr/meanptCentNbs` 下保存 POI-ref、ref-ref、Pure 和 mean-pT 的 `TProfile3D`，用于 v2-pT-rho 后处理。
- 默认输出文件为 `myAnalysisResultV2PtRho.root`。

### `calculate_c22deltapt.cpp`

- 对应 `processDataC22DeltaPt(...)`，输出与 O2 相同的 `c22DeltaPt` 目录和对象路径，也是 `c22deltapt` 四个 ROOT 后处理脚本的 AMPT 输入来源。
- charged 对象是 `c22dmeanptCharged` 和 `hMeanPtCharged`；每种 PID 粒子保存 `hMeanPt<Species>`、`c22dmeanpt<Species>RefRef`，并按模式保存 `c22dmeanpt<Species>POIRef` 或 `c22dmeanpt<Species>Pure`。
- `usePure=-1` 读取 `config.json` 中的 `c22_delta_pt_output.use_pure`，`usePure=0` 使用 PID POI-ref，`usePure=1` 使用 PID POI-POI Pure。该模式必须与 O2 的 `cfgC22DeltaPtUsePure` 含义保持一致。
- 默认输出文件为 `myAnalysisResultC22DeltaPt_pt04.root`。

### `calculate_pidptcorrelations.cpp`

- 对应 `processPidPtCorrelations(...)`，输出 `FlowContainerPidPtCorr`。
- profile 名与 O2 保持一致：`meanPtPi/Ka/Pr`，同种粒子的 `ptProductPiPi/KaKa/PrPr` 及 pair sample mean-pT，以及不同粒子的 `ptProductPiKa/PiPr/KaPr` 和相应两类粒子的 pair sample mean-pT。
- 同种粒子 pair 显式排除 self-pair；不同粒子 pair 使用两类粒子权重乘积。AMPT 权重为 `1`，但公式和 FlowContainer 的 pair weight 语义与 O2 入口一致。
- 默认输出文件为 `myAnalysisResultPidPtCorrelations.root`。

### `calculate_flowGfwOmegaXi.cpp`

- 对应 `flowGfwOmegaXi.cxx` 的 `processData(...)`，输出目录为 `flow-gfw-omega-xi`。
- 主要物理输出与 O2 `processData` 同名同轴：charged 的 `hPhi/hPhicorr/hEta/hPt/hCent/hMult/hVtxZ/hEventCount` 等基础图、`c22/c24/c32/c22Full`，以及 charged 的 `c22dpt/c24dpt/c22Fulldpt`；V0/级联粒子的 `Xic22dpt/Xic24dpt/Xic22Fulldpt/Xic24_gapdpt/Xic32dpt`、`Omegac*`、`K0sc*`、`Lambdac*` dpt profiles，以及 `InvMassXi/InvMassOmega/InvMassK0s/InvMassLambda/InvMassALambda`（含 `_all` 的 THnSparseF）。
- 从 PDG code 重建短寿命粒子：K0s 用 π+π-，Lambda/antiLambda 用 pπ-/pbarπ+，Xi/antiXi 用 π-Λ/π+antiLambda，Omega/antiOmega 用 K-Λ/K+antiLambda，并用不变质量窗口、快度窗和竞争质量排除；charged reference/POI 只要求 pT、η 和带电 PDG，不做 track-quality 筛选。
- AMPT truth 无 detector、track-quality、NUA/NUE、CCDB、local-density 等输入，因此不输出 QA、NUA 校正、MC density/locden、Jackknife 或 run-by-run detector 图；`hVtxZ/hMultTPC/hNTracksPVvsCentrality/hmultFV0AvsmultFT0A/hInteractionRate` 仅建立与 O2 同名的结构，truth 中不填充或只做最小占位。
- 默认输出文件为 `myAnalysisResultFlowGfwOmegaXi.root`。

## 参数、配置和运行

- 四个函数的共同参数依次是 input-file JSON 路径、输出 ROOT 路径、每个 JSON entry 最多读取的文件数、最多读取的 JSON entry 数、分析配置 JSON 路径；两个数量参数为 `-1` 时表示全部处理。`calculate_c22deltapt(...)` 在分析配置路径之前额外有一个 `usePure` 参数，`calculate_flowGfwOmegaXi(...)` 与 pidFlowPtCorr 三个入口的参数顺序一致。
- `config/cent_cfg.json` 列出 AMPT 输入路径前缀和文件数量；`config/config.json` 管理 flow eta subevent、独立的 mean-pT eta 范围、各粒子 pT 范围、impact-parameter 中心度转换、随机种子、Pure/POI-ref 模式以及各入口的输出轴。
- `v2_pt_rho_output`、`c22_delta_pt_output` 和 `pid_pt_correlations_output` 分别控制三个 pidFlowPtCorr 对应入口的边界策略和输出 axes。`calculate_flowGfwOmegaXi.cpp` 的输出轴按 O2 `flowGfwOmegaXi` 默认值硬编码（charged pT 38 bins、V0/级联 pT 14 bins、中心度 10 bins、mass bins 等），以保持与 O2 输出结构完全一致；修改 O2 默认 binning 时必须同步修改该脚本。修改 binning 时要保证其 O2 对应输出及 downstream 脚本仍兼容。
- 入口宏应从 `script/` 目录在 O2Physics 环境中由 ROOT 直接运行，不要编译。

```bash
cd /home/huinaibing/git_repo/amptDataAnalysis/script
root -l -b -q 'calculate_v2ptrho.cpp("../config/cent_cfg.json", "/tmp/v2ptrho.root", 1, 1, "../config/config.json")'
root -l -b -q 'calculate_c22deltapt.cpp("../config/cent_cfg.json", "/tmp/c22.root", 1, 1, -1, "../config/config.json")'
root -l -b -q 'calculate_pidptcorrelations.cpp("../config/cent_cfg.json", "/tmp/pidpt.root", 1, 1, "../config/config.json")'
root -l -b -q 'calculate_flowGfwOmegaXi.cpp("../config/cent_cfg.json", "/tmp/flowGfwOmegaXi.root", 1, 1, "../config/config.json")'
```

Codex 运行这些命令时，每条命令都必须在同一 shell 中加载 O2 环境，且 `--work-dir` 指向 `sw`：

```bash
alienv --work-dir /home/huinaibing/o2_workdir/sw \
  setenv O2Physics/latest -c bash -lc 'cd /home/huinaibing/git_repo/amptDataAnalysis/script && root -l -b -q macro.cpp'
```

不要编译 ROOT 代码；优先直接运行 ROOT 宏、脚本或命令，只有用户明确要求时才编译。

## 通用要求

- 修改前阅读相关代码和配置，遵循项目已有的风格和运行方式。
- 不要覆盖、删除或回退与当前任务无关的用户改动。
- 完成修改后进行与改动风险匹配的检查；如果无法运行某项检查，明确说明原因。
