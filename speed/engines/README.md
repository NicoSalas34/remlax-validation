# remlax engines

Time and accuracy of each remlax engine (dense, low memory, blocks, the mixed
model equations with dense, sparse and hybrid inversion) across model
families and sizes, on CPU and GPU. The protocol is [PROTOCOL.md](PROTOCOL.md).

| Folder | Content |
|---|---|
| `scripts/run_case.R` | fits one family and size with every engine |
| `scripts/families.R`, `families_params.R` | the model families |
| `scripts/analyse.R`, `analyse_params.R` | tables and figures from the results |
| `scripts/dimensions.R` | the sizes of each family |

No run is stored here. The runs of October 2026 stayed on the cluster, in
`/scratch/users/salasn/banc_moteurs_20261022`. The jobs are in
[cluster/speed_engines/](../../cluster/speed_engines/).
