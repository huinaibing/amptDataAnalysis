# Standalone O2Physics dependencies

This directory vendors only the classes used by the AMPT macros: `GFW`,
`GFWCumulant`, `GFWPowerArray`, `ProfileSubset`, and `FlowContainer`.
Their `.h`/`.cxx` files were copied from `PWGCF/GenericFramework/Core/`
in ALICE O2Physics, revision `0943ca8c368ca02e56d52d9a73e2d51f985640ce`.
Copyright notices are retained; the upstream GPL v3 license is in `COPYING`.

The only changes to the copied implementations are:

- `FlowContainer.h` includes the local `AxisSpec.h` instead of
  `Framework/HistogramSpec.h`, and uses `ClassDefInline(FlowContainer, 2)`
  to support ROOT I/O without building a dictionary library. The class name,
  persistent members, and version remain the same.
- `FlowContainer.cxx` uses ROOT's `Error()` for its one O2 `LOGF` call.

`AxisSpec.h` contains the vector/fixed-axis constructors and data members
needed by `FlowContainer::Initialize`. It is a subset of
`Framework/Core/include/Framework/HistogramSpec.h` from ALICE O2 revision
`909c5f9a70041c79be4260cba857410814613784`, with invalid fixed-axis bounds
reported by a standard exception rather than an O2 logger. It has no
dependency on O2's configurable, histogram registry, or logger classes.

`GFWRoot.h` and `FlowContainerRoot.h` include the implementation files once
per ROOT session. Use these headers in interpreted ROOT macros, without the ACLiC
`+` suffix. A ROOT 6 installation with C++17 support is required.

For downstream macros using a plain ROOT installation, include
`FlowContainerRoot.h` before reading `FlowContainer` objects. For merging,
use `script/merge_results.cpp`, which also registers the upstream `Merge`
callback with ROOT, including nominal and bootstrap profiles. The parallel
runner uses this macro automatically.
Loading these headers and the installed O2Physics versions of the same
classes together in one session is not supported; start a fresh plain ROOT
session for these standalone macros.
