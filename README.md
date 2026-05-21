# PolyhedralTemplateMatching

`PolyhedralTemplateMatching` classifies atoms using PTM and exports the reconstructed state consumed by downstream DXA-compatible tools.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- PolyhedralTemplateMatching
```

## CLI

Usage:

```bash
polyhedral-template-matching <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--crystal_structure <type>` | No | Input crystal structure: `SC`, `FCC`, `HCP`, `BCC`, `CUBIC_DIAMOND`, `HEX_DIAMOND`. | `FCC` |
| `--rmsd <float>` | No | RMSD threshold for PTM. | `0.1` |
| `--dissolve_small_clusters` | No | Mark small clusters as `OTHER` after clustering. | `false` |
| `--help` | No | Print CLI help. | |
