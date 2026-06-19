#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that extracts and clusters mitochondrial DELETION junctions from split reads.

Reads SAM (with header) on STDIN; for each primary alignment carrying an SA:Z: tag it
reconstructs the two aligned segments, and — when both segments are on the same contig
and the same strand (a deletion signature) — emits the junction:
    bp5 = end of the upstream segment  (last retained base before the deletion)
    bp3 = start of the downstream segment (first retained base after the deletion)
    SVLEN = bp3 - bp5 - 1
Junctions are clustered (breakpoints within -pad bp) and reported once per cluster with
the number of distinct supporting reads.

    EXAMPLE:
        samtools view -h -q 20 \$O.bam | sa2del.pl -chrM chrM -mtlen 16569 > \$O.sv.jun

    OUTPUT (TSV): #CHROM  bp5  bp3  SVLEN  JR  strand
~;

MAIN:
{
    my %o = (chrM => "chrM", mtlen => 16569, minsize => 50, maxsize => 0,
             pad => 10, minsupport => 2);
    GetOptions(
        "chrM=s"       => \$o{chrM},
        "mtlen=i"      => \$o{mtlen},
        "minsize=i"    => \$o{minsize},
        "maxsize=i"    => \$o{maxsize},
        "pad=i"        => \$o{pad},
        "minsupport=i" => \$o{minsupport},
        "help"         => sub { print $HELP; exit 0 },
    ) or die "ERROR: bad options\n$HELP";
    $o{maxsize} = $o{mtlen} - 1 unless ($o{maxsize});

    my @J;   # list of [up, dn, readid, strand]

    while (<>) {
        next if (/^\@/);
        my @F = split /\t/;
        next if (@F < 11);
        my $flag = $F[1];
        next unless ($flag =~ /^\d+$/);
        next if ($flag & 0x4);     # unmapped
        next if ($flag & 0x100);   # secondary
        next if ($flag & 0x800);   # supplementary (use the primary + its SA tag)
        next unless (/\tSA:Z:(\S+)/);
        my $satag = $1;

        # this (primary) segment
        my $ref    = $F[2];
        my $strand = ($flag & 0x10) ? "-" : "+";
        my ($abeg, $aend) = ref_span($F[3], $F[5]);
        next unless (defined $aend);

        # first SA segment:  rname,pos,strand,CIGAR,mapQ,NM;
        my ($sa) = split /;/, $satag;
        my @S = split /,/, $sa;
        next if (@S < 4);
        my ($sref, $spos, $sstrand, $scig) = @S[0,1,2,3];
        next unless ($sref eq $ref);          # both on the same contig
        next unless ($sstrand eq $strand);    # same orientation => deletion (not INV)
        my ($bbeg, $bend) = ref_span($spos, $scig);
        next unless (defined $bend);

        # order upstream (smaller start) / downstream
        my ($up, $dn);
        if ($abeg <= $bbeg) { $up = $aend; $dn = $bbeg; }
        else                { $up = $bend; $dn = $abeg; }
        my $svlen = $dn - $up - 1;
        next if ($svlen < $o{minsize} || $svlen > $o{maxsize});

        # distinct read identifier (mate-aware), like sam2bedSA.pl
        my $rid = $F[0];
        $rid .= "/1" if ($flag & 0x40);
        $rid .= "/2" if ($flag & 0x80);
        push @J, [$up, $dn, $rid, $strand];
    }

    # cluster by (up,dn) within pad
    @J = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @J;
    my @clusters;
    foreach my $j (@J) {
        if (@clusters && abs($j->[0] - $clusters[-1]{su}) <= $o{pad}
                      && abs($j->[1] - $clusters[-1]{sd}) <= $o{pad}) {
            push @{$clusters[-1]{pts}}, $j;
        } else {
            push @clusters, { su => $j->[0], sd => $j->[1], pts => [$j] };
        }
    }

    print "#CHROM\tbp5\tbp3\tSVLEN\tJR\tstrand\n";
    foreach my $c (@clusters) {
        my (%up, %dn, %rid, %str);
        foreach my $p (@{$c->{pts}}) {
            $up{$p->[0]}++; $dn{$p->[1]}++; $rid{$p->[2]} = 1; $str{$p->[3]}++;
        }
        my $jr = scalar keys %rid;
        next if ($jr < $o{minsupport});
        my $bp5    = mode(\%up);
        my $bp3    = mode(\%dn);
        my $strand = mode(\%str);
        my $svlen  = $bp3 - $bp5 - 1;
        next if ($svlen < $o{minsize} || $svlen > $o{maxsize});
        print join("\t", $o{chrM}, $bp5, $bp3, $svlen, $jr, $strand), "\n";
    }
    exit 0;
}

# reference span [start,end] (1-based, inclusive) consumed by a CIGAR at pos
sub ref_span {
    my ($pos, $cigar) = @_;
    return (undef, undef) unless ($pos && $pos =~ /^\d+$/ && $cigar && $cigar ne "*");
    my $len = 0;
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($n, $op) = ($1, $2);           # capture before the inner match clobbers $1
        $len += $n if ($op =~ /[MDN=X]/);  # reference-consuming ops
    }
    return ($pos, $pos + $len - 1);
}

# most frequent key (ties -> smallest numeric / lexicographic)
sub mode {
    my ($h) = @_;
    my @k = sort { $h->{$b} <=> $h->{$a} || ($a <=> $b || $a cmp $b) } keys %$h;
    return $k[0];
}
