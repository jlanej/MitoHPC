# Samplot the common deletion (del4977) for a given BAM — batteries included

A copy‑paste recipe to render a [samplot](https://github.com/ryanlayer/samplot) of the **common
mitochondrial deletion** (`del4977`, **m.8470_13447**, ~4977 bp) for any chrM BAM, using the
**MitoHPC image** — which already ships samplot (with the required matplotlib pin), `samtools`, and
the gene‑annotation track, so there is nothing to install.

> Why the image and not `pip install samplot`: a plain install pulls `matplotlib >= 3.7`, which
> regressed samplot's axes (a spurious `0..1` axis overprints the genomic axis, so the calls appear
> not to line up with the reads — samplot issues #189/#201). The image pins `matplotlib==3.6.3`, so
> the plot is correct. See `CLAUDE.md` / the `Dockerfile`.

---

## Docker

```bash
# 1. your BAM (in the current directory; adjust the -v mount otherwise) and its contig name
BAM=your.chrM.bam        # the file
MT=chrM                  # the mtDNA contig name in this BAM (see "other contig names" below)

# 2. render del4977.png from inside the MitoHPC image
docker run --rm -v "$PWD:/data" ghcr.io/jlanej/mitohpc:sv-calling bash -c "
  samtools index /data/$BAM
  samplot plot \
    -n 'del4977 common deletion (m.8470_13447)  $BAM' \
    -b /data/$BAM \
    -o /data/del4977.png \
    -c $MT -s 8470 -e 13447 -t DEL \
    -A /MitoHPC/RefSeq/genes.bed.gz
"
# -> writes del4977.png next to your BAM
```

## Singularity / Apptainer (HPC)

**Bind your data directory explicitly** with `-B` and use the bound path — robust on any site config:

```bash
BAM=your.chrM.bam
MT=chrM
singularity exec -B "$PWD:/data" docker://ghcr.io/jlanej/mitohpc:sv-calling bash -c "
  samtools index /data/$BAM
  samplot plot \
    -n 'del4977 common deletion (m.8470_13447)  $BAM' \
    -b /data/$BAM -o /data/del4977.png \
    -c $MT -s 8470 -e 13447 -t DEL \
    -A /MitoHPC/RefSeq/genes.bed.gz
"
```

> Stock Singularity/Apptainer *does* bind‑mount the current directory and set it as the working dir by
> default (`mount cwd = yes`, `--pwd`), so a bare `-b $BAM` with no `-B` often works — but that default
> is config‑dependent (frequently changed on HPC), is disabled by `--contain`/`-c` or `--no-mount cwd`,
> and can silently not apply when the data is on a scratch/network filesystem. The explicit `-B` above
> avoids all of that. (`apptainer` and `singularity` are interchangeable here.) If your data lives
> elsewhere, bind that path instead, e.g. `-B /scratch/me/run:/data`.

---

## What the command does (the "batteries")

| piece | purpose |
|-------|---------|
| `samtools index …` | samplot needs a `.bai`/`.csi`; the image's `samtools` builds it (chrM BAMs are tiny — re‑indexing is instant). |
| `-c $MT -s 8470 -e 13447 -t DEL` | draw a **deletion** between the canonical common‑deletion breakpoints **m.8470_13447**. The breakpoints actually sit anywhere inside the flanking **13 bp direct repeat** `ACCTCCCTCACCA` (8470–8482 / 13447–13459), so the bar is approximate at base resolution — that's expected. |
| `-A /MitoHPC/RefSeq/genes.bed.gz` | gene/tRNA/rRNA annotation track under the plot (already in the image). Drawn only if the track uses the same contig name as the BAM (it carries `chrM` / `rCRS` / `RSRS`). |
| `-o /data/del4977.png` | the output image. |

## How to read it

- **Coverage** (filled track): a real del4977 shows a **drop across 8470–13447** proportional to its
  heteroplasmy; a flat/continuous profile means no real deletion there.
- **Reads** (line segments): a true junction shows a **cluster of split / discordant pairs at the two
  breakpoints**. Scattered, inconsistent pairs with no coverage drop are typically artifacts.

---

## Other contig names

If your BAM does not call the mitochondrion `chrM` (e.g. `MT`, `NC_012920.1`, `rCRS`), set `MT`
accordingly. To discover it (the mtDNA contig is the one of length 16569):

```bash
docker run --rm -v "$PWD:/data" ghcr.io/jlanej/mitohpc:sv-calling \
  bash -c "samtools view -H /data/$BAM | grep LN:16569 | grep -oE 'SN:[^[:space:]]+'"
# -> e.g. SN:chrM   (use the part after 'SN:' as $MT)
```

If `$MT` is not `chrM`/`rCRS`/`RSRS`, the gene track won't match — drop the `-A …` argument (the read
and coverage panels are unaffected).

---

## Notes

- **Use the same BAM MitoHPC's caller reads.** Any chrM BAM works, but the cleanest picture comes from
  the circular‑aware, deduplicated alignment MitoHPC builds. That alignment is normally deleted after a
  run; persist it with `HP_SV_KEEPBAM=1` (writes `<prefix>.sv.bam`) and point `-b` at that.
- **Different region?** Same command, change `-s`/`-e` (and `-t DEL`/`DUP`/`INV`). For an arbitrary
  window add `-w 300` for buffer on each side.
- **Whole cohort, automatically.** A normal `mitohpc-batch-container.sh` run already produces a samplot
  gallery in `sv.report.html` for the visualizable calls (`HP_SV_PLOT`, on by default); this doc is for
  rendering a single, targeted del4977 view on demand. See `scripts/svplot.sh`.
