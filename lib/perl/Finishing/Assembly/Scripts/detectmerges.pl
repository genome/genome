#!/usr/bin/env genome-perl
#./detectmerges.pl ace_file_extension
use strict;
use warnings;
use lib "/gscuser/jschindl/mergescripts2/yfeng";
use Carp::Assert;
use Carp;
use Cwd;
use Utility;
use GSC::IO::Assembly::Ace;
use GSC::IO::Assembly::PhdDB;
use GSC::IO::Assembly::ContigTools;

if(!defined $ARGV[2])
{
	my @files = `ls  *.$ARGV[0]`;
	chomp @files;
	my ($prefix,$extension) = $files[0] =~ /(.*)\d+\.($ARGV[0])/;
	my ($outext,$num) = $extension =~ /(.+)\.(\d+)$/;
	my $outextension;
	if(defined $num)
	{
		$num++; 
		$outextension = "$outext.$num";	
	}
	else
	{
		$outextension = "$extension.2";
	}

	for(my $i = 0;$i<@files;$i++)
	{
		`bsub -o  mergelog$i -e errorlog$i detectmerges.pl $prefix$i.$extension $prefix$i.$outextension $i`;
	}
	exit;
}

my $ct = GSC::IO::Assembly::ContigTools->new;
my $ao = GSC::IO::Assembly::Ace->new(input_file => $ARGV[0],conserve_memory => 0);
my $po = GSC::IO::Assembly::PhdDB->new;

my @contig_names = @{$ao->get_contig_names};
my %super_contigs;
my @super_contig_names;
my $outfh = IO::File->new(">mergestats$ARGV[2]");
my %params = ( statsfh => $outfh, cutoffs => { hqlength => 1000, hq_percent_identity => 90 });
create_hash(@contig_names);
foreach my $super_contig_name (@super_contig_names)
{
	my $right_contig_name = $super_contigs{$super_contig_name}[0]; 
	
		
	for(my $i=1;$i<@{$super_contigs{$super_contig_name}};$i++)
	{
		my $left_contig_name = $right_contig_name;
		$right_contig_name = $super_contigs{$super_contig_name}[$i];		
		my $right_contig = $ao->get_contig($right_contig_name,1);
		my $left_contig = $ao->get_contig($left_contig_name,1);
		$ct->merge($left_contig, $right_contig,$po, %params);
		print $outfh $left_contig_name." ".$right_contig_name."\n";
		print $outfh "\n\n";
	
	}
}
#$ao->write_file(output_file => $ARGV[1]);# if ($ao->has_changed == 1);
print $outfh "Finished detecting merges\n";

#subroutines
sub create_hash
{
	my @contig_names = sort { cmp_contigs($a, $b) } @_;
	
	my $more_than_one = 0;
	foreach(@contig_names) 
	{
		if(exists $super_contigs{get_supercontig_name($_)})
		{ 
			push @{$super_contigs{get_supercontig_name($_)}}, $_;
			$more_than_one++;
			if($more_than_one == 1)
			{
				push @super_contig_names, get_supercontig_name($_);
			}
		}
		else
		{
			$super_contigs{get_supercontig_name($_)} = [$_];
			$more_than_one = 0;
		} 
	}
	return (\@super_contig_names, \%super_contigs);
}

sub get_num
{
    my ($name) = @_;
    my ($ctg_num) = $name =~ /Contig(\d+)/;	
	
    return $ctg_num;
}

sub get_supercontig_name
{
	my ($name) = @_;
	my ($ctg_num) = $name =~ /(Contig\d+)/;
	return $ctg_num;
}

sub get_ext
{
	my ($name) = @_;
	my ($ctg_ext) = $name =~ /Contig\d+\.(\d+)/;
	return $ctg_ext;
}

sub cmp_contigs
{
	my ($a, $b) = @_;
	my $num1 = get_num($a);
	my $num2 = get_num($b);
	if($num1 > $num2)
	{
		return 1;
	}elsif($num2 > $num1)
	{
		return -1;
	}
	my $ext1 = get_ext($a);
	my $ext2 = get_ext($b);
	if($ext1 > $ext2)
	{
		return 1;
	}elsif($ext2 > $ext1)
	{
		return -1;
	}
	return 0;
	
}


