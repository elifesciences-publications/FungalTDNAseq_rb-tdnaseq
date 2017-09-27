#!/usr/bin/perl -w
# Given TnSeq data from a library of transposon insertions with random barcodes, for each usable read,
# identify the barcode and the location in the genome.
#
use strict;
use Getopt::Long;
use FindBin qw($Bin);

my $minQuality = 10; # every nucleotide in a barcode must be at least this quality
my $flanking = 5; # number of nucleotides on each side that must match
my $wobbleAllowed = 2; # uncertainty in location of barcode or end of transposon, on either side of expectation
my $tmpdir = defined $ENV{TMPDIR} ? $ENV{TMPDIR} : "/tmp";
my $minIdentity = 90; # minimum %identity for mapping to genome or past-end
my $minScore = 15; # minimum score for mapping to genome or past-end
my $debug = undef;
my $trim = 0;

# Given BLAT rows (as list of lists) and hitsPastEnd hash (as reference),
# output the mapping for the read, or not (if "trumped" by hit-past-end).
# Updates hitsPastEnd for the read to score=0 if genomic hit trumps hit-past-end.
sub HandleGenomeBLAT($$);

# Global so that HandleGenomeBLAT() can use them
my $nMapped = 0;
my $nMapUnique = 0;
my $nPastEndIgnored = 0; # weak hit to past-end ignored
my $nPastEndTrumps = 0; # hit to past-end (close to) as good as hit to genome
my %nameToBarcode = ();

# Note -- should also add an option for tweaking tileSize and stepSize? 12/4 seems better than the default 11/11?

my $usage = <<END
Usage: MapTnSeq.pl [ -debug ] [ -limit maxReads ] [ -minQuality $minQuality ] [-flanking $flanking]
            [ -minIdentity $minIdentity ] [ -minScore $minScore ]
            [-tmpdir $tmpdir ]
            -genome fasta_file -model model_file -first fastq_file > output_file

    The fastq file should have phred+33 ("sanger") encoding of quality scores
    (as in MiSeq or 2012+ HiSeq). If it is named .gz, it will be gunzipped before reading.
    If it contains paired-end reads, the second read (name matching " 2:") will be ignored.

    The model file contains 1-2 lines -- the first shows what a
    typical read should look like up until the junction with the genome, e.g.
    
    nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT

    All characters in the model read must be ACGT except for the
    optional block of n\'s at the front and a block of Ns which
    represents the barcode.

    The second line in the model file is optional and contains the
    sequence "past the end" of the transposon that might arise from
    residual intact plasmid.

    This script does not handle sample multiplexing.

    The output file is tab-delimited and contains, for each usable
    read, the read name, the barcode, which scaffold the insertion
    lies in, the position of the insertion, the strand that the read
    matched, a boolean flag for if this mapping location is unique or
    not, the beginning and end of the hit to the genome in the read
    after trimming the transposon sequence, the bit score, and the
    %identity.

    minQuality specifies the minimum quality score for each character
    in the barcode.  flanking specifies the minimum number of
    nucleotides on each side that must match exactly.
END
    ;

sub FindBarcode($$$$$);
sub FindModelEnd($$$);
sub FindSubstr($$$$);
sub BLAT8($$$$$$); # BLAT to a blast8 format file

{
    my $limit = undef;
    my $fastqFile = undef;
    my $genomeFile = undef;
    my $modelFile = undef;
    my $blatcmd = -e "$Bin/blat" ? "$Bin/blat" : "blat";
    
    my $minGenomeId = 90;

    (GetOptions('debug' => \$debug, 'limit=i' => \$limit, 'minQuality=i' => \$minQuality,
                'flanking=i' => \$flanking, 'minIdentity=i' => \$minIdentity, 'minScore=i' => \$minScore, 'trim=i' => \$trim,
                'tmpdir=s' => \$tmpdir,
                'blat=s' => \$blatcmd,
                'genome=s' => \$genomeFile, 'model=s' => \$modelFile, 'first=s' => \$fastqFile)
     && @ARGV==0) || die $usage;
    die $usage unless defined $modelFile && defined $fastqFile;

    die "Cannot read $genomeFile" unless -r $genomeFile;
    die "Cannot read $fastqFile" unless -r $fastqFile;
    die "Not a directory: $tmpdir" unless -d $tmpdir;
    die "minScore must be at least 10" if $minScore < 10;
    die "Invalid minimum identity" if $minIdentity < 50 || $minIdentity > 100;

    open(MODEL, "<", $modelFile) || die "Cannot read $modelFile";
    my $model = <MODEL>;
    $model =~ s/[\r\n]+$//;

    die "Invalid model: $model" unless $model =~ m/^n*[ACGT]+N+[ACGT]+$/;
    my $pastEnd = <MODEL>;
    if (defined $pastEnd) {
	chomp $pastEnd;
	die "Invalid past-end sequence: $pastEnd" unless $pastEnd =~ m/^[ACGT]+$/;
}	
    close(MODEL) || die "Error reading $modelFile";

    my $barcodeStart = index($model, "N");
    my $barcodeEnd = rindex($model, "N");
    my $barcodeLen = $barcodeEnd-$barcodeStart+1;
    print STDERR "Parsed model $modelFile\n";
    print STDERR "Barcodes of length $barcodeLen, expected transposon region of " . length($model);
    print STDERR ", no pastEnd" if !defined $pastEnd;
    print STDERR "\n";

    $model = substr($model,0,-$trim);

    my $pipe = 0;
    if ($fastqFile =~ m/[.]gz$/) {
        $pipe = 1;
        open(FASTQ, '-|', 'zcat', $fastqFile) || die "Cannot run zcat on $fastqFile";
    } elsif ($fastqFile =~ m/[.zip]$/) {
        $pipe = 1;
        open(FASTQ, "7za -so e $fastqFile | zcat |") || die "Cannot run 7za and zcat on $fastqFile";
    } else {
        open(FASTQ, "<", $fastqFile) || die "Cannot read $fastqFile";
    }

    my %nameToHits = (); # list of scaffold, position, strand, match score
    my $nLong = 0;

    my $rand = rand();
    my $tmpFna = $tmpdir . "/MapTnSeq_" . $$ . "_$rand.fna";
    open(TMPFNA, ">", $tmpFna) || die "Cannot write to $tmpFna";

    my $nReads = 0;
    my $nTryToMap = 0;

    # Find barcodes and end of transposon and write remaining sequence to TMPFNA
    while(!defined $limit || $nReads < $limit) {
        my $name = <FASTQ>;
        (last, next) unless $name;
        chomp $name;
        
        my $seq = <FASTQ>;
        chomp $seq;
        <FASTQ>;
        my $quality = <FASTQ>;
        chomp $quality;

        next if $name =~ m/^\S+ 2:/; # ignore second side of paired-end reads
        $nReads++;

        # short sequences are unmappable
        next unless length($seq) >= length($model) + $minScore;
        $nLong++;

        my ($barcode,$obsStart) = FindBarcode($seq,$quality,$model,$barcodeStart,$barcodeEnd);
        next unless defined $barcode;

        $name =~ s/ .*$//;
        die "Duplicate read name: $name" if exists $nameToBarcode{$name};
        $nameToBarcode{$name} = $barcode;

        my $transposonEnd = FindModelEnd($seq,$model,$obsStart - $barcodeStart);
        if (defined $transposonEnd && length($seq) >= $transposonEnd + $minScore) {
            print STDERR "Try to map $name\n" if $debug;
            my $inGenome = substr($seq, $transposonEnd+1+$trim);
            print TMPFNA ">$name\n$inGenome\n";
            $nTryToMap++;
        }
    }

    # Prematurely closing the pipe gives an error
    close(FASTQ) || ($pipe && defined $limit) || die "Error reading from $fastqFile: $!";
    close(TMPFNA) || die "Error writing to $tmpFna";
    print STDERR "Read $nReads reads\n";

    my %hitsPastEnd = (); # read to score of match to past-end sequence
    my %hitsPastEndDetails = (); # read to score of match to past-end sequence

    if (defined $pastEnd) {
	# Map to past-end-of transposon
	my $endFna = $tmpdir . "/MapTnSeq_end" . $$ . "_" . "_$rand.fna";
	open(END, ">", $endFna) || die "Cannot write to $endFna";
	print END ">pastend\n$pastEnd\n";
	close(END) || die "Error writing to $endFna";
	my $blat8 = BLAT8($tmpFna, $endFna, $tmpdir, $blatcmd, $minScore, $minIdentity);
	print STDERR "Parsing past-end hits to $blat8\n" if defined $debug;
	open(BLAT, "<", $blat8) || die "Cannot read $blat8";
	while(<BLAT>) {
	    chomp;
	    my @F = split /\t/, $_;
	    my ($query, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $eval, $score) = @F;
	    my $oldscore = 0;
	    if (exists $hitsPastEnd{$query}) {
		my @oldHitsPastEnd = split(":",$hitsPastEnd{$query});
		$oldscore = $oldHitsPastEnd[0];
	    }
	    if ($oldscore < $score) {
		$hitsPastEnd{$query} = join(":",($score,$sBeg,$sEnd,$identity,$qBeg,$qEnd));
            }
	}
	close(BLAT) || die "Error reading $blat8";
	unlink($blat8) unless defined $debug;
    }

    # Map to the genome
    my $blat8 = BLAT8($tmpFna, $genomeFile, $tmpdir, $blatcmd, $minScore, $minIdentity);
    print STDERR "Parsing $blat8\n" if defined $debug;
    open(BLAT, "<", $blat8) || die "Cannot read $blat8";
    my @lines = ();
    while(<BLAT>) {
        chomp;
        my @F = split /\t/, $_;
        my ($query, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $eval, $score) = @F;
        if (@lines == 0 || $query eq $lines[0][0]) {
            push @lines, \@F;
        } else {
            HandleGenomeBLAT(\@lines, \%hitsPastEnd);
            @lines = \@F;
        }
    }
    HandleGenomeBLAT(\@lines, \%hitsPastEnd);

    close(BLAT) || die "Error reading $blat8";
    unlink($blat8) unless defined $debug;

    unlink($tmpFna) unless defined $debug;

    # print out hits-past-end
    while (my ($read,$blastResult) = each %hitsPastEnd) {

	my ($score,$sBeg,$sEnd,$identity,$qBeg,$qEnd) = split(":",$blastResult);

        next unless $score > 0; # score=0 means it was ignored above
        print join("\t", $read, $nameToBarcode{$read}, "pastEnd", $sBeg, ($sBeg < $sEnd ? "+" : "-"), 1, $qBeg, $qEnd, $score, $identity)."\n";

    }

    print STDERR "Reads processed $nReads Long-enough $nLong Barcodes found " . scalar(keys %nameToBarcode) . "\n";
    print STDERR "Mapping attempted for $nTryToMap Mapped $nMapped Uniquely $nMapUnique\n";
    print STDERR "Hits past end of transposon: " . (scalar(keys %hitsPastEnd) - $nPastEndIgnored) .
        " plus $nPastEndIgnored weak/ambiguous; trumped hit to genome $nPastEndTrumps times\n";
    print STDERR sprintf("Proportions: Long-enough %.3f Barcode %.3f Attempted %.3f Mapped %.3f Past-end %.3f\n",
			 $nLong/$nReads,
			 scalar(keys %nameToBarcode)/$nReads,
			 $nTryToMap/$nReads,
			 $nMapped/$nReads,
			 (scalar(keys %hitsPastEnd) - $nPastEndIgnored)/$nReads)
	if $nReads > 0;
}

sub FindSubstr($$$$) {
    my ($subseq, $seq, $expAt, $wobble) = @_;
    my $len = length($seq);
    for (my $i = $expAt - $wobble; $i <= $expAt+$wobble; $i++) {
        return($i) if $i >= 0 && $i < $len && substr($seq, $i, length($subseq)) eq $subseq;
    }
    return undef;
}

sub FindBarcode($$$$$) {
    my ($seq,$quality,$model,$expStart,$expEnd) = @_;
    my $pre = substr($model, $expStart - $flanking, $flanking);
    my $post = substr($model, $expEnd + 1, $flanking);
    my $preLoc = FindSubstr($pre, $seq, $expStart - $flanking, $wobbleAllowed);
    return undef unless defined $preLoc;
    my $postLoc = FindSubstr($post, $seq, $expEnd+1, $wobbleAllowed);
    return undef unless defined $postLoc;

    # positions of 1st and last nucleotides of barcode
    my $start = $preLoc + $flanking;
    my $end = $postLoc-1;

    my $barcodeLength = ($end-$start);
    
    unless (abs(($expEnd-$expStart) - ($end-$start)) < 4) {
        print STDERR "Wrong barcode length: $start $end not $expStart $expEnd in $seq\n" if defined $debug;
        return undef;
    }
    my $barcode = substr($seq, $start, $end-$start+1);

    if ($minQuality > 0) {
        # the sanger code for !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
        my $barqual = substr($quality, $start, $end-$start+1);
        # quality score is encoded as 33+
        my @scores = map { $_ - 33 } unpack("%C"x$barcodeLength, $barqual);
        foreach my $score (@scores) {
            die "Invalid score $score from barcode $barcode quality $barqual" if $score < 0 || $score > 100;
            if ($score < $minQuality) {
                print STDERR "Low quality $score for barcode $barcode in $seq\n" if defined $debug;
                return undef;
            }
        }
    }
    return($barcode,$start);
}

sub FindModelEnd($$$) {
    my ($seq,$model,$expOff) = @_;
    my $expEnd = length($model) + $expOff;
    my $at = FindSubstr(substr($model,length($model)-$flanking,$flanking), $seq, $expEnd-$flanking, $wobbleAllowed);
    if (!defined $at) {
        print STDERR "No end sequence " . substr($model,length($model)-$flanking,$flanking) . " near position " . ($expEnd-$flanking) . " of\n$seq\n" if defined $debug;
        return undef;
    }
    return($at + $flanking - 1); # last position of transposon
}

# returns the name of the file
sub BLAT8($$$$$$) {
    my ($queriesFile,$dbFile,$tmpdir,$blatcmd,$minScore,$minIdentity) = @_;

    my $blat8File = "$tmpdir/MapTnSeq.$$.".rand() . ".psl";
    # prevent writing to stdout by BLAT
    open(OLD_STDOUT, ">&STDOUT");
    print OLD_STDOUT ""; # prevent warning
    open(STDOUT, ">", "/dev/null");
    system($blatcmd, "-out=blast8", "-t=dna", "-q=dna", "-minScore=$minScore", "-minIdentity=$minIdentity", "-maxIntron=0", "-noTrimA", $dbFile, $queriesFile, $blat8File) == 0
        || die "Cannot run $blatcmd: $!";
    close(STDOUT);
    open(STDOUT, ">&OLD_STDOUT") || die "Cannot restore stdout: $!";
    return($blat8File);
}

sub HandleGenomeBLAT($$) {

    my ($rows, $hitsPastEnd) = @_;
    return if scalar(@$rows) == 0;
    my $read = $rows->[0][0];
    die if !defined $read || !exists $nameToBarcode{$read};

    my @hits = ();

    # indexes within besthits entries
    my ($SCAFFOLD,$POSITION,$STRAND,$SCORE,$QBEG,$QEND) = (0,1,2,3,4,5);
    my @besthits = ();

    foreach my $row (@$rows) {
        my ($read2, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $eval, $score) = @$row;
        die unless $read2 eq $read;
        if (scalar(@besthits) == 0 || $score >= $besthits[0][$SCORE] - 5) {
            # convert from 0-based to 1-based position, and note that sBeg always < sEnd so flip if stranded
            push @besthits,  [ $subject, $sBeg, ($sBeg < $sEnd ? "+" : "-"), $score,
                               $identity, $qBeg, $qEnd ];
            print STDERR "Hit for $read:$qBeg:$qEnd to $subject:$sBeg:$sEnd score $score identity $identity\n"
                if $debug;
        }
    }

    die if @besthits == 0;
    my ($scaffold, $position, $strand, $score, $identity, $qBeg, $qEnd) = @{ $besthits[0] };

    # and output a mapping row (or none)
    if (exists $hitsPastEnd->{$read}) {

	my ($HPEscore,$HPEsBeg,$HPEsEnd,$HPEidentity,$HPEqBeg,$HPEqEnd) = split(":",$hitsPastEnd->{$read});
        if ($HPEscore >= $score - 5) {
            $nPastEndTrumps++;
            print STDERR "Past end trumps for $read\n" if $debug;
            return;
        } else {
            $nPastEndIgnored++; # ignore weak hit to past-end sequence
            print STDERR "Ignoring hit past end of transposon for $read\n" if $debug;
            $hitsPastEnd->{$read} = join(":",(0,$HPEsBeg,$HPEsEnd,$HPEidentity,$HPEqBeg,$HPEqEnd)); # so we do not print out a line for it later on
        }
    }
    # else if no trumping
    print join("\t", $read, $nameToBarcode{$read}, $scaffold, $position, $strand, @besthits == 1 ? 1 : 0,
               $qBeg, $qEnd, $score, $identity)."\n";
    $nMapUnique++ if @besthits == 1;
    $nMapped++;
    return;
}
