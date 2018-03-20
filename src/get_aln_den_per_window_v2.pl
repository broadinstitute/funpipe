#!/usr/bin/env perl
# Get Alignment Cover in defined Sliding Window Size
#	Run Command:
#	get_aln_den_per_window.pl

my ($mPileup, $winSize) = @ARGV;
my $winEnd = $winSize - 1;
my $currSliWin = 1;
my $currContig = 1;
my $depthSum = 0;
my $loop = 0;
my $avgWin = 0;
my $currWin;
if ($mPileup =~ /gz$/) {
	open (MPILEUP, "zcat $mPileup |") || die "can't open ALN file\n";
} else {
	open (MPILEUP, "<$mPileup") || die "can't open ALN file\n";
}

while (<MPILEUP>)
{
	my $line = $_;
	chomp $line;
	my ($contig, $position, $depth) = split "\t", $line;

	if($currContig ne $contig)
	{
		$currWin = ($currSliWin*$winSize) - $winSize;
		print "$currContig\t$currWin\t$currSliWin\t$depthSum\n";

		# set the current contig window end and rowcount
		$currContig = $contig;
		$winEnd = $winSize -1;
		$currSliWin = 1;
		$depthSum = 0;
		$avgWin = 0;
	}
	else
	{
		while($loop==0)
		{
			if ($position < $winEnd)
			{
				$depthSum = $depthSum + $depth;
				$loop = 1;
			}
			else
			{
				$avgWin = ($depthSum/5000);
				$currWin = ($currSliWin*$winSize)-$winSize;
				print "$currContig\t$currWin\t$currSliWin\t$depthSum\n";
				$depthSum = 0;
				$currSliWin = $currSliWin + 1;
				$winEnd = $winEnd + $winSize;
				$avgWin = 0;
			}
		}
		$loop = 0;
	}
}

$currWin = ($currSliWin*$winSize)-$winSize;
print "$currContig\t$currWin\t$currSliWin\t$depthSum\n";
