FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

ENV HP_HDIR=/MitoHPC/
ENV HP_SDIR=/MitoHPC/scripts/
ENV HP_BDIR=/MitoHPC/bin/
ENV HP_ADIR="bams"
ENV HP_ODIR="out"
ENV HP_IN="in.txt"

ENV PATH="$HP_SDIR:$HP_BDIR:$PATH"

###########################################

RUN apt-get -y update
# minimap2: used by test/sv/gen_bams.sh to (re)generate the SV mock BAMs in-image.
# The SV caller itself needs only samtools/bedtools/perl (installed below); the
# committed mock BAMs let test/sv/run_test.sh run without minimap2.
RUN apt-get install -y wget tar nano curl git minimap2

###########################################
# HEAVY, SLOW-CHANGING LAYER (~25 min: whole-genome reference download + tool compiles + jars).
# Built from ONLY the install machinery + the committed RefSeq/ — deliberately NOT the SV caller,
# tests, or docs — so editing callsv.py / callSV.sh / the tests / docs is a CACHE HIT here and
# skips the whole download+compile (locally and via the GHA layer cache in CI).
#
# Keep this COPY list in sync with what the install RUN below invokes:
#   init.sh                     sourced for HP_* vars (its ls2in.pl line is guarded by
#                               `[ -d bams/ ]`, which is absent at build time, so it never runs)
#   install_sysprerequisites.sh apt packages + `pip install pysam` (so checkInstall's import works)
#   install_prerequisites.sh    tools/jars/references; calls circFasta.sh + rotateFasta.sh via PATH
#   circFasta.sh, rotateFasta.sh  self-contained (samtools/bedtools/inline-perl); only run if a
#                               circularized/rotated reference is missing
#   checkInstall.sh             verifies tools/jars/refs/pysam
# RefSeq/ is copied so install_prerequisites.sh's `[ ! -s rCRS.dict ]`-style guards SKIP rebuilding
# the references (they are committed); only the uncommitted whole-genome hs38DH.fa is downloaded.
COPY scripts/init.sh scripts/install_sysprerequisites.sh scripts/install_prerequisites.sh \
     scripts/circFasta.sh scripts/rotateFasta.sh scripts/checkInstall.sh /MitoHPC/scripts/
COPY RefSeq/ /MitoHPC/RefSeq/
RUN \
  chmod a+x $HP_SDIR/*.* && \
  . $HP_SDIR/init.sh && \
  $HP_SDIR/install_sysprerequisites.sh && \
                                                                     $HP_SDIR/install_prerequisites.sh && \
  export HP_MT=rCRS && export HP_MTC=rCRSC && export HP_MTR=rCRSR && $HP_SDIR/install_prerequisites.sh && \
  export HP_MT=RSRS && export HP_MTC=RSRSC && export HP_MTR=RSRSR && $HP_SDIR/install_prerequisites.sh && \
  HP_MT=RSRS $HP_SDIR/install_prerequisites.sh && \
  $HP_SDIR/checkInstall.sh && \
  rm -fr /MitoHPC/prerequisites/

# samplot — structural-variant visualization for the SV module's optional HP_SV_PLOT step. TWO
# version constraints both matter for correct, stable plots:
#   * matplotlib==3.6.3 — matplotlib>=3.7 REGRESSED samplot's axes: a spurious 0..1 normalized axis
#     is overprinted on the real genomic x-axis and the coverage/insert-size y-axes, so breakpoints
#     become unreadable and calls appear not to line up with the reads (samplot issues #189/#201; the
#     maintainers' fix is matplotlib<3.7, 3.6.3 known-good). The bug is glaring on narrow (small-
#     deletion) windows and subtle on wide ones. numpy<2 mirrors the known-good set and keeps cp310
#     wheels (no compiler). This — NOT the samplot commit — was the actual cause of the bad plots; the
#     earlier unpinned `pip install matplotlib` + `--no-deps` let mpl>=3.7 in and bypassed samplot's
#     own matplotlib<3.7 guard.
#   * samplot @2929e4a — the v1.3.0 release code (identical plot rendering) PLUS the Python-3.11/numpy
#     forward-compat fixes (v1.3.0 itself crashes on py>=3.11: random.sample over dict_keys). PyPI
#     ships a broken 0.0.1, so install from git. --no-deps so samplot can't pull a matplotlib over the
#     pins above; jinja2 is its remaining runtime dep (pysam is already present for callsv.py).
# The build-time assert fails the image if matplotlib ever resolves to >=3.7. Cache-stable.
RUN python3 -m pip install --no-cache-dir 'matplotlib==3.6.3' 'numpy<2' jinja2 && \
    python3 -m pip install --no-cache-dir --no-deps \
      'samplot @ git+https://github.com/ryanlayer/samplot.git@2929e4a' && \
    python3 -c 'import matplotlib as m; v=tuple(map(int,m.__version__.split(".")[:2])); assert v<(3,7), "matplotlib "+m.__version__+" reintroduces the samplot axis-label bug (#189/#201)"' && \
    samplot plot --help >/dev/null

###########################################
# FAST LAYER: bring in the rest of the repo (SV caller, tests, docs). A code-only edit re-runs
# only from here, reusing the cached install layer above. Re-chmod (the full scripts/ is present
# now) and drop the frozen example fixtures. The SV caller is smoke-tested against the committed
# mock data by CI (.github/workflows/docker-publish.yml runs `run_test.sh` in the built image), so
# it is not re-run inside the build.
COPY . /MitoHPC/
RUN chmod a+x $HP_SDIR/*.* && rm -fr /MitoHPC/examples*
