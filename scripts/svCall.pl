#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that corroborates clustered deletion junctions with read depth, estimates
heteroplasmy two independent ways, applies false-positive flags, and emits a VCF body
(STDOUT) plus a flat table (-tab).

A deletion is PASS only when it has both >= -minjr junction reads AND a real coverage
drop (medInside/medFlank <= -drop) AND enough flanking depth, and is not origin-wrapped.

    INPUT (STDIN, from sa2del.pl): #CHROM bp5 bp3 SVLEN JR strand
    EXAMPLE:
        sa2del.pl ... | svCall.pl -sample S -ref chrM.fa -cvg S.dp -mtlen 16569 \\
            -hp HP.bed.gz -numt NUMT.vcf.gz -dloop DLOOP.bed.gz -tab S.sv.tab > body.vcf

Heteroplasmy:
    AFJ = JR/(JR+SR), SR = max(round(mean(depth[bp5],depth[bp3])) - JR, 0)   (junction fraction)
    AFC = 1 - medInside/medFlank                                            (coverage ratio)
    AFDIFF = |AFJ-AFC|   (QC: large => amplification bias or DUP-as-DEL)
~;

MAIN:
{
    my %o = (sample => "SAMPLE", mtlen => 16569, flank => 200, minjr => 3,
             drop => 0.9, mindepth => 0, originpad => 20, pad => 25,
             rep5a => 8470, rep5b => 8482, rep3a => 13447, rep3b => 13459);
    my ($cvg, $ref, $hp, $numt, $dloop, $tab);
    GetOptions(
        "sample=s"   => \$o{sample}, "mtlen=i"   => \$o{mtlen},
        "flank=i"    => \$o{flank},  "minjr=i"   => \$o{minjr},
        "drop=f"     => \$o{drop},   "mindepth=i"=> \$o{mindepth},
        "originpad=i"=> \$o{originpad}, "pad=i"  => \$o{pad},
        "cvg=s" => \$cvg, "ref=s" => \$ref, "hp=s" => \$hp,
        "numt=s" => \$numt, "dloop=s" => \$dloop, "tab=s" => \$tab,
        "help" => sub { print $HELP; exit 0 },
    ) or die "ERROR: bad options\n$HELP";
    die "ERROR: -cvg required\n" unless ($cvg);

    # depth array (1-based)
    my @dep = (0) x ($o{mtlen} + 1);
    open(my $C, "<", $cvg) or die "ERROR: cannot read $cvg: $!";
    while (<$C>) { my @f = split; next unless (@f >= 3 && $f[1] =~ /^\d+$/);
                   $dep[$f[1]] = $f[2] if ($f[1] <= $o{mtlen}); }
    close($C);

    # reference sequence (for the VCF REF base)
    my $seq = "";
    if ($ref && -s $ref) {
        open(my $R, "<", $ref) or die;
        while (<$R>) { next if (/^>/); chomp; $seq .= $_; }
        close($R);
    }

    my @hpiv    = load_bed($hp);
    my @dliv    = load_bed($dloop);
    my %numtpos = load_vcf_pos($numt);

    open(my $T, ">", $tab) or die "ERROR: cannot write $tab: $!" if ($tab);
    print $T "#sample\tchrom\tbp5\tbp3\tsvlen\tJR\tSR\tAFJ\tAFC\tAFDIFF\tCVGR\tFLANKDP\tFILTER\tflags\n" if ($tab);

    while (<STDIN>) {
        next if (/^#/);
        chomp;
        my ($chrom, $bp5, $bp3, $svlen, $jr, $strand) = split /\t/;
        next unless (defined $jr && $bp5 =~ /^\d+$/);

        my $medIn   = median_range(\@dep, $bp5 + 1, $bp3 - 1, $o{mtlen});
        my $medFl   = median_flank(\@dep, $bp5, $bp3, $o{flank}, $o{mtlen});
        my $ratio   = $medFl > 0 ? $medIn / $medFl : 1;
        my $span    = int(($dep[wrap($bp5,$o{mtlen})] + $dep[wrap($bp3,$o{mtlen})]) / 2 + 0.5);
        my $sr      = $span - $jr; $sr = 0 if ($sr < 0);
        my $afj     = ($jr + $sr) > 0 ? $jr / ($jr + $sr) : 0;
        my $afc     = 1 - $ratio; $afc = 0 if ($afc < 0); $afc = 1 if ($afc > 1);
        my $afdiff  = abs($afj - $afc);

        # flags
        my @flags;
        my $repeat = (near($bp5, $o{rep5a}, $o{rep5b}, $o{pad}) || near($bp3, $o{rep3a}, $o{rep3b}, $o{pad}) ||
                      near($bp5, $o{rep3a}, $o{rep3b}, $o{pad}) || near($bp3, $o{rep5a}, $o{rep5b}, $o{pad}));
        my $wrap   = ($bp5 <= $o{originpad} || $bp5 >= $o{mtlen} - $o{originpad} ||
                      $bp3 <= $o{originpad} || $bp3 >= $o{mtlen} - $o{originpad});
        my $inhp   = in_bed(\@hpiv, $bp5) || in_bed(\@hpiv, $bp3);
        my $indl   = in_bed(\@dliv, $bp5) || in_bed(\@dliv, $bp3);
        my $innumt = (exists $numtpos{$bp5} || exists $numtpos{$bp3});
        push @flags, "REPEAT" if ($repeat);
        push @flags, "WRAP"   if ($wrap);
        push @flags, "HP"     if ($inhp);
        push @flags, "DLOOP"  if ($indl);
        push @flags, "NUMT"   if ($innumt);

        # PASS / FILTER
        my @fil;
        push @fil, "lowJR"       if ($jr < $o{minjr});
        push @fil, "no_cvg_drop" if ($ratio > $o{drop});
        push @fil, "WRAP"        if ($wrap);
        push @fil, "lowDP"       if ($o{mindepth} && $medFl < $o{mindepth});
        my $filter = @fil ? join(";", @fil) : "PASS";

        my $refbase = ($bp5 >= 1 && $bp5 <= length($seq)) ? substr($seq, $bp5 - 1, 1) : "N";
        $refbase = "N" if ($refbase eq "");
        my $end  = $bp3 - 1;
        my $info = sprintf("SM=%s;SVTYPE=DEL;END=%d;SVLEN=%d;JR=%d;SR=%d;AFJ=%.3f;AFC=%.3f;AFDIFF=%.3f;CVGR=%.3f",
                           $o{sample}, $end, -$svlen, $jr, $sr, $afj, $afc, $afdiff, $ratio);
        $info .= ";$_" foreach (@flags);

        printf("%s\t%d\t.\t%s\t<DEL>\t.\t%s\t%s\tGT:DP:AF\t0/1:%d:%.3f\n",
               $chrom, $bp5, $refbase, $filter, $info, int($medFl + 0.5), $afj);

        if ($tab) {
            printf $T "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%s\t%s\n",
                $o{sample}, $chrom, $bp5, $bp3, $svlen, $jr, $sr, $afj, $afc, $afdiff,
                $ratio, int($medFl + 0.5), $filter, (@flags ? join(",", @flags) : ".");
        }
    }
    close($T) if ($tab);
    exit 0;
}

sub wrap { my ($p, $m) = @_; $p = (($p - 1) % $m) + 1; return $p; }

sub median_range {
    my ($dep, $a, $b, $m) = @_;
    return 0 if ($b < $a);
    my @v;
    for (my $p = $a; $p <= $b; $p++) { push @v, $dep->[wrap($p, $m)]; }
    return median(\@v);
}

sub median_flank {
    my ($dep, $bp5, $bp3, $flank, $m) = @_;
    my @v;
    for (my $p = $bp5 - $flank + 1; $p <= $bp5; $p++) { push @v, $dep->[wrap($p, $m)]; }
    for (my $p = $bp3; $p <= $bp3 + $flank - 1; $p++) { push @v, $dep->[wrap($p, $m)]; }
    return median(\@v);
}

sub median {
    my ($v) = @_;
    return 0 unless (@$v);
    my @s = sort { $a <=> $b } @$v;
    my $n = @s;
    return ($n % 2) ? $s[int($n/2)] : ($s[$n/2 - 1] + $s[$n/2]) / 2;
}

sub near { my ($p, $a, $b, $pad) = @_; return ($p >= $a - $pad && $p <= $b + $pad); }

# load BED.gz (chrom start end ; 0-based half-open) -> list of [start1,end1] 1-based inclusive
sub load_bed {
    my ($f) = @_;
    my @iv;
    return @iv unless ($f && -s $f);
    open(my $H, "gzip -dc '$f' |") or return @iv;
    while (<$H>) { next if (/^#/); my @x = split; next unless (@x >= 3 && $x[1] =~ /^\d+$/);
                   push @iv, [$x[1] + 1, $x[2]]; }
    close($H);
    return @iv;
}

sub in_bed {
    my ($iv, $p) = @_;
    foreach my $r (@$iv) { return 1 if ($p >= $r->[0] && $p <= $r->[1]); }
    return 0;
}

# load VCF.gz positions (handles any chrom)
sub load_vcf_pos {
    my ($f) = @_;
    my %h;
    return %h unless ($f && -s $f);
    open(my $H, "gzip -dc '$f' |") or return %h;
    while (<$H>) { next if (/^#/); my @x = split; next unless (@x >= 2 && $x[1] =~ /^\d+$/); $h{$x[1]} = 1; }
    close($H);
    return %h;
}
