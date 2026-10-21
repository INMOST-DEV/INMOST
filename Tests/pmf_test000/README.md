# PMF MPI-IO regression

`pmf_test000` exercises distributed hexahedra with shared entities, one ghost
layer and a real tag on nodes and cells. It compares the complete output of
strategies 0 and 1, loads MPI-IO output using both readers, checks counts and tag
values, and overwrites the same path with smaller meshes. Finally, it removes
the last payload byte and requires `BadFile` on **every** rank. Only one rank
receives a short read in that case; error handling must remain collective.

Build the `pmf_test000` target with `USE_MESH`, `USE_MPI`, `USE_MPI_FILE` and
`COMPILE_TESTS` enabled, then run:

```sh
ctest --test-dir <build> -R pmf_test000 --output-on-failure
```

## Cluster evidence, 2026-10-21

The failure was observed on the
[Sechenov HPC cluster](https://sechenov-hpc.readthedocs.io/ru/latest/?badge=latest)
using INMOST commit `2896b43e6f326d6610d58e38f59386b48ea337b8`. The workload
repeatedly saved and immediately loaded a distributed hexahedral mesh with
12 cells per axis per rank, shared entities, one ghost layer and a real tag on
nodes and cells. All four comparison runs used 32 ranks on two nodes, with
16 ranks per node and NFSv4 storage. Each run requested 50 iterations.

The observed results were:

- Strategy 1 (MPI-IO), new directories: `BadFile` loading iteration 5.
- Strategy 0 (rank-zero file I/O with MPI transfers): all 50 iterations
  succeeded on all 32 ranks.
- Strategy 1, existing directory: `BadFile` loading iteration 2.
- Strategy 1, one-second delay before Load: `BadFile` at iteration 1.

Iteration indices above are zero-based. Save returned on all ranks before the
failed Load calls.

Inspection of the PMF file from the first failed run showed an incomplete
32-part file:

- Its prefix is `10 01 08`; the outer header uses little-endian 8-byte integers.
- The part count is 32, so the outer header occupies 267 bytes.
- Part 0 and part 31 declare 794874 bytes each. Parts 1 through 30 declare
  865542 bytes each. Expected total size: **27556275 bytes**.
- Actual size: **13778271 bytes**, exactly `267 + 794874 + 15 * 865542`.
- Parts 0 through 15 have their mesh headers, original rank identifiers and
  final `09` markers at the expected offsets. Parts 16 through 31 are absent.
- The first missing byte is the start of rank 16's part, precisely the node
  boundary. Error probes in all three failed runs report this same file size.

This establishes missing payloads in the inspected file, not just an invalid
prefix or an immediate-read timing issue. It does not establish whether the
second node wrote a different/unlinked file identity or whether MPI/NFS lost or
exposed stale file size/data. No cluster MPI implementation trace is available.
The fixed MPI message tags and `USE_MPI_IRECV_ORDER_RULE` do not govern the
strategy-1 file writes; their changes in 2896b43e are in the strategy-0 branch.
The strategy-1 write sequence already existed before that commit.

## Change and limits

The writer now opens once, sets the exact final size collectively and separates
resizing from data writes with sync/barrier/sync. This removes the unnecessary
create/delete/recreate cycle and avoids extending the file independently from
multiple nodes. It also truncates stale tails when overwriting a larger file.
Writes are flushed on every rank before any rank closes its handle.

The reader checks actual payload transfer counts, closes the file collectively
and propagates a short-read failure to all ranks before reconstruction. An
MPI-IO EOF can return `MPI_SUCCESS` with fewer bytes; previously the reader
treated the unread buffer tail as file contents and ranks could diverge.

MPI references: [resizing](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node355.htm),
[consistency](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node383.htm),
[sync examples](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node393.htm).

Local filesystem tests alone cannot reproduce or prove elimination of the
two-node NFS failure. Additional sync operations have a potential performance
cost; quantify it with comparable successful runs on the cluster storage.

Local validation: Debug, MPICH, `USE_MPI_IRECV_ORDER_RULE=ON`; mesh enabled,
solver/partitioner/autodiff/nonlinear modules disabled. All six selected CTest
cases (`pmf_test000`, `mpi_tag_test000`, `mpi_window_test000`, each with 2 and 4
ranks) passed. A separate run of the workload described above also completed
50 iterations with 4 ranks, 12 cells per axis per rank, strategy 1, new
per-iteration directories and no delay.

Cluster validation after the change used the same workload on the Sechenov HPC
cluster: 32 ranks on two nodes, 12 cells per axis per rank, strategy 1, NFSv4,
50 iterations, new per-iteration directories and no delay. All 1600
Save calls and 1600 Load calls returned, and all 32 ranks finalized successfully.
This confirms that the previously failing scenario passes with the change.
The successful run does not isolate
which part of the writer change resolves the underlying MPI/NFS behavior and
does not provide a before/after performance comparison.
