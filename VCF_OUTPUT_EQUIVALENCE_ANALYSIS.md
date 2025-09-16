# MitoHPC VCF Output Equivalence Analysis

## Question
Will `mitohpc-parallel.sh` produce the same sets of merged VCFs as `mitohpc.sh`? Such as `mutect2.mutect2.10.merge.vcf`

## Answer: YES

The merged VCF files will be **identical** between both approaches. Here's the detailed analysis:

## Processing Workflow Comparison

### Sequential Processing (`mitohpc.sh`)
1. `run.sh` generates processing commands
2. Commands executed sequentially via `bash`
3. Each sample processed in order from `$HP_IN`
4. `getSummary.sh` called after all samples complete

### Parallel Processing (`mitohpc-parallel.sh`)
1. `run.sh` generates processing commands → `run.all.sh`
2. Extract `filter.sh` commands → `filter.commands.sh`
3. Run `filter.sh` commands **in parallel** (any order)
4. Run `getSummary.sh` after all samples complete

## Key Finding: Output Order vs Processing Order

The critical insight is that **processing order does not affect output order** because:

### 1. Sample Order Determination
Both scripts use the same sample order from `$HP_IN` file:
```bash
# From getSummary.sh line 35:
awk '{print $3}' $HP_IN | sed "s|$|.$S.00.vcf|" | xargs cat
```

### 2. VCF Concatenation Process
- Individual sample VCF files are concatenated in `$HP_IN` order (NOT processing order)
- The concatenated file is sorted by genomic coordinates: `bedtools sort -header`
- Merged VCF is generated from the sorted concatenated file

### 3. Merge VCF Generation Pipeline
```bash
# From snpCount.sh line 69:
cat $HP_ODIR/$S.$T.concat.vcf | concat2merge.pl -in $HP_IN | tee $HP_ODIR/$S.$T.merge.vcf | vcf2sitesOnly.pl > $HP_ODIR/$S.$T.merge.sitesOnly.vcf
```

## Why Outputs Are Identical

1. **Same Input Files**: Both approaches process the exact same BAM/CRAM files
2. **Same Sample Order**: VCF concatenation follows `$HP_IN` order regardless of processing order
3. **Same Sorting**: `bedtools sort -header` ensures genomic coordinate sorting
4. **Same Processing Logic**: Identical filtering, merging, and output generation scripts

## Performance vs Consistency Trade-off

- **mitohpc.sh**: Sequential processing, slower but preserves processing order
- **mitohpc-parallel.sh**: Parallel processing, faster but processes samples in any order

**Result**: Same output files, different processing time.

## Verification Approach

To verify this analysis, compare the merged VCF files from both approaches:

```bash
# Run both approaches on the same dataset
mitohpc.sh
mitohpc-parallel.sh 4

# Compare merged VCF files
diff out/mutect2.mutect2.10.merge.vcf out_parallel/mutect2.mutect2.10.merge.vcf
# Should show no differences
```

## Conclusion

**Yes, `mitohpc-parallel.sh` will produce the same sets of merged VCFs as `mitohpc.sh`**. The parallelization affects only the processing speed, not the output consistency. The merged VCF files like `mutect2.mutect2.10.merge.vcf` will be identical between both approaches.