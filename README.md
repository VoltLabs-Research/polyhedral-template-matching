# Polyhedral Template Matching

Runs PTM, exports the per-atom structure type, and generates the cluster-graph artifacts consumed by OpenDXA.

## Install

```bash
vpm install @voltlabs/polyhedral-template-matching
```

## CLI

```bash
polyhedral-template-matching <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--crystal_structure <type>` | no | `FCC` | Input crystal structure: `SC`, `FCC`, `HCP`, `BCC`, `CUBIC_DIAMOND`, `HEX_DIAMOND`. |
| `--rmsd <float>` | no | `0.1` | RMSD threshold for PTM (min `0`). |
| `--dissolve_small_clusters` | no | `false` | Mark small clusters as `OTHER` after clustering. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_atoms.parquet` | Structure Identification | AtomisticExporter → glb |
| `{output_base}_atoms.parquet` | Structure Counts Chart | ChartExporter → chart-png |
| `{output_base}_ptm_analysis.parquet` | PTM Analysis | — |
| `{output_base}_clusters.table` | Clusters Table | — |
| `{output_base}_cluster_transitions.table` | Clusters Transitions | — |
| `{output_base}_atoms.parquet` | Per Atom Properties | — |
| `{output_base}_neighbor_lattice.parquet` | Neighbor Lattice | — |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/polyhedral-template-matching
