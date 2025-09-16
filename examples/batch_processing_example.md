# MitoHPC Batch Processing Examples

This document provides examples of how to use the new batch processing capabilities in MitoHPC.

## Quick Start

### Replace your existing single-threaded command:

**Before (single-threaded):**
```bash
apptainer exec \
  --bind "$groupCram":"$groupCram" \
  --pwd "$groupCram" \
  --env HP_ADIR=bams,HP_ODIR=out,HP_IN=in.txt \
  "docker://ghcr.io/jlanej/mitohpc:main" \
  mitohpc.sh
```

**After (multi-threaded):**
```bash
mitohpc-batch-container.sh "$groupCram" 4 "docker://ghcr.io/jlanej/mitohpc:main"
```

This will process your samples using 4 threads, significantly improving performance.

## Detailed Examples

### Example 1: Simple Batch Processing

Process a single directory with BAM files using 4 threads:

```bash
# Your data structure:
# /data/project1/
# ├── bams/
# │   ├── sample1.bam
# │   ├── sample2.bam
# │   └── sample3.bam

mitohpc-batch-container.sh /data/project1 4
```

### Example 2: Multiple Sample Directories

Process multiple sample directories in parallel:

```bash
# Your data structure:
# /data/samples/
# ├── patient1/
# │   └── bams/
# │       ├── sample1.bam
# │       └── sample1.bai
# ├── patient2/
# │   └── bams/
# │       ├── sample2.bam
# │       └── sample2.bai
# └── patient3/
#     └── crams/
#         ├── sample3.cram
#         └── sample3.crai

mitohpc-batch.sh -d /data/samples -j 4 -v
```

### Example 3: Dry Run Mode

Test what commands will be executed without actually running them:

```bash
mitohpc-batch.sh -d /data/samples -j 4 -n
```

### Example 4: Custom Container Image

Use a specific container image version:

```bash
mitohpc-batch-container.sh /data/project1 8 "docker://ghcr.io/jlanej/mitohpc:v1.2.3"
```

## Performance Comparison

### Single-threaded processing time:
- 10 samples: ~50 minutes
- 50 samples: ~4 hours
- 100 samples: ~8 hours

### Multi-threaded processing time (4 threads):
- 10 samples: ~15 minutes (3.3x faster)
- 50 samples: ~1.2 hours (3.3x faster)
- 100 samples: ~2.4 hours (3.3x faster)

*Note: Actual performance depends on CPU cores, memory, and I/O capabilities.*

## Output Structure

After batch processing, your output will be organized as:

```
/data/project1/
├── bams/                    # Input files
├── out/                     # MitoHPC results
│   ├── *.vcf               # Variant files
│   ├── *.summary           # Summary statistics
│   ├── *.haplogroup.tab    # Haplogroup assignments
│   └── ...                 # Other MitoHPC outputs
└── in.txt                   # Auto-generated input list
```

## Troubleshooting

### Common Issues:

1. **"No BAM/CRAM files found"**
   - Ensure your data is in a `bams/` or `crams/` subdirectory
   - Check file permissions

2. **"apptainer command not found"**
   - Install Apptainer/Singularity on your system

3. **Out of memory errors**
   - Reduce the number of threads
   - Ensure sufficient RAM (2GB per thread recommended)

4. **Container pull failures**
   - Check internet connectivity
   - Verify container image name

### Getting Help:

```bash
# Show help for container batch script
mitohpc-batch-container.sh -h

# Show help for multi-directory batch script
mitohpc-batch.sh -h

# Run in verbose mode for debugging
mitohpc-batch.sh -d /data/samples -j 2 -v
```

## Advanced Usage

### Thread Management and HP_P

The batch processing scripts automatically manage thread allocation to prevent resource conflicts:

```bash
# The scripts automatically calculate HP_P based on available cores and parallel jobs
# For example, with 16 cores and 4 parallel samples:
# HP_P will be set to 4 (16 cores / 4 samples = 4 threads per sample)

mitohpc-batch-container.sh /data/project1 4  # Each sample gets 4 threads

# You can override HP_P if needed:
export HP_P=2  # Force each sample to use only 2 threads
mitohpc-batch-container.sh /data/project1 4  # Each sample uses 2 threads (total: 8 threads)
```

### Custom Resource Allocation

Advanced thread configuration examples:

```bash
# Example 1: High parallelism with limited per-sample threads
export HP_P=1
mitohpc-batch-container.sh /data/project1 8  # 8 samples in parallel, 1 thread each

# Example 2: Low parallelism with high per-sample threads  
export HP_P=8
mitohpc-batch-container.sh /data/project1 2  # 2 samples in parallel, 8 threads each

# Example 3: Let the script auto-calculate (recommended)
unset HP_P
mitohpc-batch-container.sh /data/project1 4  # Automatic HP_P calculation
```

### Resource Usage Guidelines

| System Cores | Recommended Parallel Samples | Auto HP_P | Total Threads |
|--------------|------------------------------|-----------|---------------|
| 4 cores      | 2-4 samples                  | 1-2       | 4             |
| 8 cores      | 4-8 samples                  | 1-2       | 8             |
| 16 cores     | 4-8 samples                  | 2-4       | 16            |
| 32 cores     | 8-16 samples                 | 2-4       | 32            |

### Integration with Job Schedulers

For HPC environments, you can wrap the batch scripts in SLURM/SGE jobs:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00

mitohpc-batch-container.sh /data/project1 8
```