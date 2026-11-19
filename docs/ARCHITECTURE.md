<!-- Managed by MetaViro-DRANK GitHub setup -->
# Repository architecture

| Directory | Purpose |
| --- | --- |
| `bin/` | Human-facing executable shell wrappers |
| `scripts/` | Actual validated Python analysis modules and orchestrators |
| `docs/` | Manuals, design and validation record |
| `examples/` | Instructions and small distributable example configurations only |
| `results/` | Generated outputs and local tests (ignored by Git) |

The existing root Python entry points are preserved for backward compatibility. Do **not** move validated modules into `bin/` just to reorganize the repository: module imports, `__file__`-relative locations, shell invocations, and provenance expectations may depend on current paths. Future structural refactoring should be a separate, tested change.

The Python programs inside `scripts/` remain the source of truth. Do not copy them into `bin/`; this would create divergent implementations.
