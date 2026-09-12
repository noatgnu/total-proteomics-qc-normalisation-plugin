# Total Proteomics QC & Normalisation


## Installation

**[⬇️ Click here to install in Cauldron](http://localhost:50060/install?repo=https%3A%2F%2Fgithub.com%2Fnoatgnu%2Ftotal-proteomics-qc-normalisation-plugin)** _(requires Cauldron to be running)_

> **Repository**: `https://github.com/noatgnu/total-proteomics-qc-normalisation-plugin`

**Manual installation:**

1. Open Cauldron
2. Go to **Plugins** → **Install from Repository**
3. Paste: `https://github.com/noatgnu/total-proteomics-qc-normalisation-plugin`
4. Click **Install**

**ID**: `total-proteomics-qc-normalisation`  
**Version**: 1.0.0  
**Category**: preprocessing  
**Author**: CauldronGO Team

## Description

QC filtering and PRONE normalisation method comparison for DIA-NN total proteomics data

## Runtime

- **Environments**: `r`

- **Entrypoint**: `total_proteomics_qc_normalisation.R`

## Inputs

| Name | Label | Type | Required | Default | Visibility |
|------|-------|------|----------|---------|------------|
| `pg_matrix_file` | Protein Group Matrix File | file | Yes | - | Always visible |
| `stats_file` | Stats File | file | No | - | Always visible |
| `annotation_file` | Sample Annotation File | file | Yes | - | Always visible |
| `min_unique_peptides` | Minimum Proteotypic Peptides | number (min: 0, step: 1) | No | 2 | Always visible |
| `contaminant_column` | Contaminant Column | text | No | Contaminant | Always visible |
| `cc_mapped_column` | Curated Category Column | text | No | - | Always visible |

### Input Details

#### Protein Group Matrix File (`pg_matrix_file`)

DIA-NN report.pg_matrix.tsv (protein group intensities, one column per raw file)


#### Stats File (`stats_file`)

DIA-NN report.stats.tsv (per-run QC metrics). Optional; omit to skip run-level QC plots.


#### Sample Annotation File (`annotation_file`)

Cauldron sample annotation file. Sample is matched against pg_matrix column names by exact string, falling back to basename.

- **Table Editor**: Enabled with 3 columns
  - **Columns**:
    - `Sample`: Sample (required)
      - Sample identifier matching a pg_matrix column (full raw-file path or basename)
    - `Condition`: Condition (required)
      - Experimental group label
    - `BioReplicate`: BioReplicate
      - Biological replicate identifier (optional)

#### Minimum Proteotypic Peptides (`min_unique_peptides`)

Minimum N.Proteotypic.Sequences required to keep a protein group


#### Contaminant Column (`contaminant_column`)

Column flagging contaminants with '+'; matching rows are dropped. Skipped automatically if pg_matrix has no such column.


#### Curated Category Column (`cc_mapped_column`)

Optional column whose value equal to the primary gene name flags a 'CC.mapped' category (annotation only, never filtered on). Leave empty to skip.


## Outputs

| Name | File | Type | Format | Description |
|------|------|------|--------|-------------|
| `sample_annotation` | `sample_annotation.tsv` | data | tsv | Resolved sample-to-group mapping used for the run |
| `log2_matrix` | `log2_matrix.tsv` | data | tsv | Filtered protein groups as a log2 intensity matrix, before normalisation |
| `normalized_median` | `normalized_Median.tsv` | data | tsv | PRONE Median-normalised matrix |
| `normalized_quantile` | `normalized_Quantile.tsv` | data | tsv | PRONE Quantile-normalised matrix |
| `normalized_vsn` | `normalized_VSN.tsv` | data | tsv | PRONE VSN-normalised matrix |
| `normalized_normicsvsn` | `normalized_NormicsVSN.tsv` | data | tsv | PRONE NormicsVSN-normalised matrix (skipped automatically if PRONE cannot fit it on this dataset) |
| `cv_comparison_plot` | `normalisation/CV_comparison_all_methods.png` | plot | png | Per-group CV across normalisation methods; the deciding plot for choosing a method for the differential-expression plugin |
| `qc_run_summary_plot` | `QC/QC_run_summary.png` | plot | png | Combined per-run QC panel (proteins, precursors, missed cleavages identified). Present only when a stats file was given. |

## Requirements

- **R Version**: >=4.0

### R Dependencies (External File)

Dependencies are defined in: `r-packages.txt`

- `SummarizedExperiment`
- `PRONE`
- `tidyverse`
- `ggplot2`
- `RColorBrewer`
- `scales`
- `svglite`
- `patchwork`

> **Note**: When you create a custom environment for this plugin, these dependencies will be automatically installed.

## Example Data

This plugin includes example data for testing:

```yaml
  pg_matrix_file: examples/pg_matrix.tsv
  stats_file: examples/stats.tsv
  annotation_file: examples/annotation.txt
  min_unique_peptides: 2
```

Load example data by clicking the **Load Example** button in the UI.

## Usage

### Via UI

1. Navigate to **preprocessing** → **Total Proteomics QC & Normalisation**
2. Fill in the required inputs
3. Click **Run Analysis**

### Via Plugin System

```typescript
const jobId = await pluginService.executePlugin('total-proteomics-qc-normalisation', {
  // Add parameters here
});
```
