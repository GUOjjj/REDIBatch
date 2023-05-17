#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# To get uniquely Mapped reads from HISAT2 BAM

my ($bam,$outdir,$suffix);
my $samtools;

GetOptions(
	"bam=s" => \$bam,
	"outdir=s" => \$outdir,
	"suffix=s" => \$suffix,
	"samtools=s" => \$samtools,
);

$samtools ||="samtools";
$suffix ||="sort.rmdup.bam";

`mkdir -p $outdir` unless (-e $outdir);

die unless (-e $samtools);

my $name=(split /\//,$bam)[-1];
$name=~s/\.$suffix$//;
my $outbam="$outdir/$name.sort.rmdup.unique.bam";
open OT,"| $samtools view -b -S ->$outbam" or die $!;
open IN,"$samtools view -h $bam|" or die $!;
while (<IN>){
	chomp;
	if (/^@/){
		print OT "$_\n";
		next;
	}
	print OT "$_\n" if (/NH:i:1$/);
}
close IN;
close OT;

