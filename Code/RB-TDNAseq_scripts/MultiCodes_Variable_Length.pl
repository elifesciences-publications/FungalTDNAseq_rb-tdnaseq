#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $minQuality = 10;

my $usage = <<END
Usage: MultiCodes.pl [ -debug ] [ -dntag ] [ -limit maxReads ] [ -minQuality $minQuality ]
          [ -index name | -primers PrimerIndexTable ]
          [ -preseq CAGCGTACG -postseq AGAGACCTC -nPreExpected 9 ]
 -out out_prefix < fastq
    PrimerIndexTable should be tab-delimited with multiplex name and a primer like nACGACG
    The fastq file should be fastq with phred+33 ("sanger") encoding of quality scores
    (as in MiSeq or 2012+ HiSeq)

This script analyzes multiplexed random barcode sequences, such as
nATCACGAG GTCGACCTGCAGCGTACG 20N AGAGACCTCGTGGACATCAGATC
for primer nATCACGAG (nLeading=1, indexseq=ATCACGAG), where the 20 Ns are the barcode

or, if -dntag is specified,
nATCACGAG CGGTGTCGGTCTCGTAG 20N CGATGAATTCGAGCTCGTT

or, if the index is specified, there is no multiplexing, and the read is
nnnnn GTCGACCTGCAGCGTACG 20N AGAGACC (where the leading ns are random and are ignored)

or, if -preseq -postseq -nPreExpected are all used, it expects a read of the form
   nPreExpected characters (any, although this includes the multiplexing tag unless -index is used)
   preseq
   barcode
   postseq

For a barcode to be counted, the preseq (9 nt upstream of the barcode)
much match exactly; the postseq (9 nt downstream of the barcode) much
also be correct unless the sequence is too short (in which 4 or 0 are
checked) or minQuality = 0 (in which case the post-sequence is not
checked and there is no guarantee that the barcode is the correct
length).
END
    ;

sub Variants($); # return all 1-nt variants for a sequence
sub FindPrefix($$); # sequence, list of index prefixes => index of prefix or -1
sub FindBarcode($$$); # sequence, quality, first post-prefix index => (barcode, spacing) or undef
sub sum(@); # sum of a list of values

my $preseq = undef; # before the barcode
my $postseq = undef; # after (for 50 nt reads this needs to be shortened to AGA)
my $nPreExpected = 0; # nt between prefix and preseq
my $debug = 0;
my $dntag = 0;
my $iname = undef;

{
    my ($indexfile,$out,$nLimit);
    (GetOptions('primers=s' => \$indexfile,
		'out=s' => \$out,
		'index=s' => \$iname,
                'minQuality=i' => \$minQuality,
                'limit=i' => \$nLimit,
                'nPreExpected=i' => \$nPreExpected,
		'preseq=s' => \$preseq,
		'postseq=s' => \$postseq,
                'dntag' => \$dntag,
                'debug' => \$debug)
     && defined $out)
        || die $usage;
    die $usage unless (defined $indexfile xor defined $iname);

    if (defined $preseq) {
	die "Missing -postseq: $usage" unless defined $postseq;
	die "Missing -nPreExpected: $usage" unless defined $nPreExpected;
    } else {
	if ($dntag) {
	    die "-index with -dntag not supported" if defined $iname;
	    $preseq = "GTCTCGTAG";
	    $nPreExpected = 8;
	    $postseq = "CGATGAATT";
	} else {
	    $preseq = "CAGCGTACG";
	    $nPreExpected = defined $iname ? 14 : 9;
	    $postseq = "AGAGACCTC";
	}
    }
    my $nReads = 0;
    my $nMulti = 0; # number with prefix identified
    my %nOff = map {$_ => 0} (18..22); # only spacings of 20 are considered; count totals
    my %codes = (); # barcode => number of times barcode is seen, with an array of one entry per multiplex tag
    my @prefix = (); # list of [ nLeading, indexseq (no leading Ns), name ]

    my @prefixNames = ();
    if (defined $indexfile) {
	open(INDEX, "<", $indexfile) || die "Error reading $indexfile";
	while (my $line = <INDEX>) {
	    chomp $line;
	    my ($name,$index) = split /\t/, $line;
	    next if $name =~ m/name/i;
	    my $nLeading;
	    my $indexseq;
	    if ($name eq "none" && $index eq "") {
		$nLeading = 0;
		$indexseq = "";
	    } else {
		die "Invalid index sequence $index" unless $index =~ m/^([nN]*)([ACGT]+)$/;
		$nLeading = length($1);
		$indexseq = $2;
	    }
	    die $line if !defined $nLeading || !defined $indexseq || !defined $name;
	    push @prefix, [ $nLeading, $indexseq, $name ];
	}
	close(INDEX) || die "Error reading $indexfile";
	print STDERR "Read " . scalar(@prefix) . " indices from $indexfile\n";
	@prefixNames = map { $_->[2] } @prefix;
    }
    if (defined $iname) {
	@prefix = ( [0, "", $iname] );
	@prefixNames = ( $iname );
    }

    while(my $header = <STDIN>) {
        chomp $header;
        my $seq = <STDIN>;
        chomp $seq;
        $nReads++;
        <STDIN>;
        my $quality = <STDIN>;

        last if defined $nLimit && $nReads >= $nLimit;

        my $iPrefix = defined $iname ? 0 : FindPrefix($seq, \@prefix);
        next if $iPrefix < 0;
        $nMulti++;
        my ($nLeading, $indexseq, $prefixName) = @{ $prefix[$iPrefix] };
        my ($barcode,$off) = FindBarcode($seq, $quality, $nLeading+length($indexseq));
        next unless defined $barcode;
        $nOff{$off}++;
        next unless ($off > 16 && $off < 22);
        if (!exists $codes{$barcode}) {
            $codes{$barcode} = [];
        }
        $codes{$barcode}[$iPrefix]++;
    }

    my $nOffTot = sum(values %nOff);
    my $nUniq = scalar(keys %codes);
    print STDERR "Reads $nReads Multiplexed $nMulti Usable(20) $nOff{20} fraction " . ($nOff{20}/$nReads) . " unique codes $nUniq \n" if $nReads > 0;
    foreach my $off (17..22) {
        print STDERR sprintf("Off\t%d\t%d\t%.3f\n", $off, $nOff{$off}, $nOff{$off}/$nOffTot) if $nOffTot > 0;
    }

    my %nPerCount = (); # number of codes with that count
    open(CODES, ">", "$out.codes.tmp") || die "Cannot write to $out.codes.tmp";
    print CODES join("\t","barcode", @prefixNames)."\n";
    while (my ($code,$count) = each %codes) {
        my @counts = ();
        my $nPrefix = scalar(@prefix);
        for (my $i=0; $i < $nPrefix; $i++) {
            push @counts, $count->[$i] || 0;
        }
        print CODES join("\t",$code,@counts)."\n";
        my $nAll = sum(@counts);
        $nPerCount{$nAll}++;
    }
    close(CODES) || die "Error writing to $out.codes.tmp";
    rename("$out.codes.tmp","$out.codes") || die "Cannot rename to $out.codes";
    print STDERR "Wrote " . scalar(keys %codes) . " unique barcodes to $out.codes\n";

    open(COUNTS, ">", "$out.counts") || die "Cannot write to $out.counts";
    print COUNTS join("\t", "Count", "nCodes", "Frac")."\n";
    foreach my $count (sort {$a<=>$b} keys %nPerCount) {
        print COUNTS join("\t",$count,$nPerCount{$count},$nPerCount{$count}/$nUniq)."\n";
    }
    close(COUNTS) || die "Error writing to $out.counts";
    print STDERR "Wrote number of barcodes seen a certain number of times to $out.counts";
    print STDERR "; nOnce = $nPerCount{1}" if (exists $nPerCount{1});
    print STDERR "\n";

    # and look for off-by-1 cases, unless minQuality is off
    if ($minQuality != 0) {
        open(CLOSE, ">", "$out.close") || die "Cannot write to $out.close";
        print CLOSE join("\t",qw{code1 count1 code2 count2})."\n";
        my $nCases = 0;

        while (my ($code,$count) = each %codes) {
            my @variants = Variants($code);
            foreach my $variant (@variants) {
                if (($code cmp $variant) > 0 && exists $codes{$variant}) {
                    print CLOSE join("\t",$code,sum(@$count),$variant,sum(@{$codes{$variant}}))."\n";
                    $nCases++;
                }
            }
        }
        close(CLOSE) || die "Error writing to $out.close";
        print STDERR "Wrote $nCases off-by-1 pairs to $out.close\n";
    }
}

sub FindPrefix($$) {
    my ($seq, $indexes) = @_;
    my @matches = ();
    for (my $i = 0; $i < scalar(@$indexes); $i++) {
        my ($nLeading,$indexseq,$name) = @{ $indexes->[$i] };
        if (substr($seq, $nLeading, length($indexseq)) eq $indexseq) {
            push @matches, $i;
        }
    }
    return $matches[0] if @matches == 1;
    print STDERR "No prefix for $seq\n" if $debug;
    return -1;
}

sub FindBarcode($$$) {
    my ($seq,$quality,$offset) = @_;
    my $seq2 = substr($seq, $offset);
    my $quality2 = substr($quality, $offset);

    my $prepos = index($seq2, $preseq);
    unless($prepos >= 0 && $prepos >= $nPreExpected-2 && $prepos <= $nPreExpected+2) {
        print STDERR "seq2 $seq2 has invalid index-of-preseq $prepos\n" if $debug;
        return undef;
    }

    my $postseqUsed = $postseq;
    $postseqUsed = "" if $minQuality == 0;
    if (length($seq2) < $prepos + length($preseq) + 20 + length($postseqUsed)) {
        $postseqUsed = substr($postseqUsed, 0, 4);
        print STDERR "Using postseq $postseqUsed\n" if $debug;

        if (length($seq2) < $prepos + length($preseq) + 20) {
            print STDERR "Giving up, too short\n" if $debug;
            return undef;
        }
        if (length($seq2) < $prepos + length($preseq) + 20 + length($postseqUsed)) {
            print STDERR "Ignoring postseq, too short\n" if $debug;
            $postseqUsed = "";
        }
    }

    my $foundEnd = -1;
    if ($postseqUsed eq "") {
        $foundEnd = 20;
    } else {
        # need to check 20 first in case end of barcode matches post-seq. (AGAG issue)
        foreach my $off (20,19,21,18,17,22) {
            if (substr($seq2, $prepos + length($preseq) + $off, length($postseqUsed)) eq $postseqUsed) {
                $foundEnd = $off;
                last;
            }
        }
    }
    if ($foundEnd == -1) {
        print STDERR "seq2 $seq2 barcode postseq not found\n" if $debug;
        return undef;
    }

    my $barcode = substr($seq2, $prepos+length($preseq), $foundEnd);
    unless(length($barcode) < 22 && length($barcode) > 16 && $barcode =~ m/^[A-Z]+$/) {
        # note barcodes with ambiguous sequences are ignored
        # (we might want to allow one N in the future)
        print STDERR "seq2 $seq2 has invalid barcode sequence $barcode\n" if $debug;
        return undef;
    }
    
    if ($minQuality > 0) {
        # the sanger code for !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
        my $barqual = substr($quality2, $prepos+length($preseq), $foundEnd);
        # quality score is encoded as 33+
        my @scores = map { $_ - 33 } unpack("%C"x$foundEnd, $barqual);
        foreach my $score (@scores) {
            die "Invalid score $score from barcode $barcode quality $barqual" if $score < 0 || $score > 100;
            if ($score < $minQuality) {
                print STDERR "Low quality $score for barcode $barcode in $seq2\n" if $debug;
                return undef;
            }
        }
    }
    return($barcode,$foundEnd);
}

sub Variants($) {
    my ($baseseq) = @_;
    my @out = ();
    $baseseq = uc($baseseq);
    my $len = length($baseseq);
    foreach my $i (0..($len-1)) {
        my $pre = substr($baseseq,0,$i);
        my $char = substr($baseseq,$i,1);
        my $post = substr($baseseq,$i+1);
        next unless $char eq "A" || $char eq "C" || $char eq "G" || $char eq "T";
        foreach my $newchar (qw{A C G T}) {
            push @out, $pre . $newchar . $post unless $newchar eq $char;
        }
    }
    return(@out);
}

sub sum( @ ) {
    my @x = @_;
    my $sum = 0;
    foreach (@x) { $sum += $_ if defined $_; }
    return $sum;
}
