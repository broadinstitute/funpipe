#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Spec::Functions 'catfile';
use Data::Dumper;

##### subroutines #####

sub help{
   print STDERR qq{
   Usage: dep_per_win.pl
      --output_dir|o       Output directory
      --prefix|p           Output prefix
      --window             Sliding window size
      --mpileup|m          Pileup
      --faidx              Fasta index file
      --h|help\n
   };
   exit 2;
}

sub output_tsv{
  # output
  my ($prefix, $out_dir, $chrDep, $winSize) = @_;
  my $output = catfile($out_dir, $prefix.".depth.tsv");
  open (OUT, "+>$output") || die "can't open out file\n";
  print OUT join("\t", qw(chr start0 end0 id depth fc)), "\n";
  foreach my $currContig (sort keys %$chrDep) {
    foreach my $currSliWin (sort {$a <=> $b} keys $$chrDep{$currContig}) {
      my $currWin = get_curr_win($currSliWin, $winSize);
      print OUT sprintf("%s\t%d\t%d\t%s\t%.2f\t%.2f\n", $currContig, $currWin,
                        $currWin + $winSize, $currContig.'_'.$currSliWin,
                        $$chrDep{$currContig}{$currSliWin}{dep},
                        $$chrDep{$currContig}{$currSliWin}{fc});
    }
  }
  close (OUT) or die;
}

sub depth_fold_change {
  # Normalize depth by overall average depth.
  my ($winDep, $aveDep) = @_;
  foreach my $contig (keys %$winDep) {
    foreach my $winIdx (keys %{$$winDep{$contig}}) {
      if ($aveDep == 0) {
        $$winDep{$contig}{$winIdx}{fc} = 'inf';
      } else {
        $$winDep{$contig}{$winIdx}{fc} =
          $$winDep{$contig}{$winIdx}{dep}/$aveDep;
      }
    }
  }
}

sub average_depth{
  # About: calculate mean depth per each chromosome
  # Usage: mean_dep_per_chr($chrlen, $chrDep)
  # Args:  $chrlen: pointer to a hash table contains contig length of each chr
  #        $chrDep: pointer to a hash table contains overall depth of coverage of
  #                 each chr
  # Return:%chrcov: a hash table contains depth of coverage for each chr.

  my ($chrDep, $chrlen) = @_;
  my ($dep, $len) = (0, 0);   # Total depth and total length
  foreach my $contig (keys %$chrDep) {
    $dep += $$chrDep{$contig};
    $len += $$chrlen{$contig};
  }
  return ($dep/$len);
}

sub load_faidx {
  # 	load fasta file index
	my ($file) = @_;
	my %chrlen;
	open(IN, $file) or die;
	while(<IN>){
		chomp;
		my ($chr, $len, $offset, $linebases, $linewidth) = split(/\t/, $_);
		$chrlen{$chr} = $len;
	}
	close(IN) or die;
	return(\%chrlen)
}

sub update_chr_dep{
  # update chromosome depth
	my ($chrDep, $contig, $depth) = @_;
	if (!defined($$chrDep{$contig})) {
		$$chrDep{$contig} = $depth;
	} else {
		$$chrDep{$contig} += $depth;
	}
}

sub get_curr_win {
  #  get current window from window size
  my ($currSliWin, $winSize) = @_;
  return $currSliWin * $winSize;
}

##### main script #####
my $nOpts = scalar(@ARGV);

# default parameter
my $winSize = 1000;
my $prefix = 'dep_per_win';
my $out_dir = '.';

&GetOptions(
   "help|h" => \my $help,
	 "output_dir|o=s" => \$out_dir,
	 "prefix|p=s" => \$prefix,
   "window=s" => \$winSize,
	 "mpileup|m=s" => \my $mPileup,
	 "faidx=s" => \my $faidx
);
&help() if $help || $nOpts == 0;

# initalize
my $winEnd = $winSize - 1;  # end position of window
my $currContig = 'NA';         # current contig
my $depthSum = 0;           # sum of depth within a window
my $loop = 0;               # loop
my $avgWin = 0;             # average window size
my $currWin;                # beginning position of current window
my $currSliWin = 0;         # 1-based index of sliding window size

my $chrDep = {};            # reference to chr depth hash
my $winDep;                 # average depth of each sliding window

if ($mPileup =~ /gz$/) {
	open (MPILEUP, "zcat $mPileup |") || die "can't open ALN file\n";
} else {
	open (MPILEUP, "<$mPileup") || die "can't open ALN file\n";
}
while (my $line = <MPILEUP>) {
	chomp $line;
	my ($contig, $position, $depth) = split "\t", $line;
	update_chr_dep($chrDep, $contig, $depth);
	if ($currContig ne $contig) {
		$currWin = get_curr_win($currSliWin, $winSize);
		# set the current contig window end and rowcount
		$currContig = $contig;
		$winEnd = $winSize -1;
		$currSliWin = 0;           # 0 based sliding window
		$depthSum = 0;
	} else {
		while ($loop==0) {
			if ($position < $winEnd) {
				$depthSum += $depth;
				$loop = 1;
			} else {
				$currWin = get_curr_win($currSliWin, $winSize);
        $$winDep{$currContig}{$currSliWin}{dep} = $depthSum/$winSize;
        # initialize
        $depthSum = 0;
        $currSliWin += 1;
        $winEnd += $winSize;
      }
    }
    $loop = 0;
  }
}
close (MPILEUP) or die;
$$winDep{$currContig}{$currSliWin}{dep} = $depthSum/$winSize;

#
my $chrlen = load_faidx($faidx);
my $aveDep = average_depth($chrDep, $chrlen);
depth_fold_change($winDep, $aveDep);
output_tsv($prefix, $out_dir, $winDep, $winSize);

__END__
=head1 NAME: dep_per_win.pl

=head1 SYNOPSIS

   ./dep_per_win.pl --help

=head1 DESCRIPTION

  Get Alignment Cover in defined Sliding Window Size

=head1 AUTHOR

  This script is modified from Jose Munoz's code by Xiao Li.
  For bugs/questions, contact xiaoli@broadinstitute.org.

=head1 LOGS

  03/26/18 16:39:03
    -  Initial implementation of code.

=cut
