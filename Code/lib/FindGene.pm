
package FindGene;
require Exporter;
use strict;
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw( LocationToGene );

# Given scaffold, pos, and hash of scaffold to list of sorted genes,
# returns the locusId and the fraction through the gene it is in (as a list of 2 elements)
#
# If the location overlaps multiple genes or no genes, returns locusId = "", f = "".
#
# Each gene should be a hash that contains begin, end, strand, and locusId
#
# This code does not support orfs that wrap around the origin.
sub LocationToGene($$$) {
    my ($scaffold, $pos, $sortedGenes) = @_;
    return ("","") if $scaffold eq "pastEnd";

    my $genelist = $sortedGenes->{$scaffold};
    return ("","") if !defined $genelist;

    # binary search
    # at all times, either the true index is between lo and hi, or there is no hit
    my $nGenes = scalar(@$genelist);
    my $lo = 0;
    my $hi = $nGenes-1;
    for (my $nRound = 0; $nRound < 100000; $nRound++) {
	my $mid = int(($lo+$hi)/2);
	my $iBegin = $genelist->[$mid]{begin};
	my $iEnd = $genelist->[$mid]{end};

	if ($pos < $iBegin) {
	    return ("","") if $mid == $lo;
	    $hi = $mid-1;
	} elsif ($pos > $iEnd) {
	    return ("","") if $mid == $hi;
	    $lo = $mid+1;
	} else {
	    # does the previous or next gene also overlap this position?
	    return ("","") if ($mid > 0 && $genelist->[$mid-1]{begin} <= $pos && $pos <= $genelist->[$mid-1]{end});
	    return ("","") if ($mid < $nGenes-1 && $genelist->[$mid+1]{begin} <= $pos && $pos <= $genelist->[$mid+1]{end});
	    my $f = $iBegin == $iEnd ? 0 : ($pos - $iBegin)/($iEnd-$iBegin);
	    my $strand = $genelist->[$mid]{strand};
	    # insertions near N terminus of gene should have f near 0 regardless of strand
	    $f = 1.0-$f if $strand eq "-";
	    return($genelist->[$mid]{locusId}, $f);
	}
    }
    die "Unreachable";
}

1;
