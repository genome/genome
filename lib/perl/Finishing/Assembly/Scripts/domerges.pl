#!/usr/bin/env genome-perl
# usage: domerges.pl ace_file_extension mergelist > mergescript
# to execute, simply run the mergescript
use IO::File;
use strict;
use warnings;
my $fh = IO::File->new($ARGV[1]);
my @files = `ls *scaffold*.$ARGV[0]`;
chomp @files;
my $mod = scalar @files;
my ($prefix) = $files[0] =~ /(.*)\d+\.$ARGV[0]/;
sub get_num
{
    my ($name) = @_;
    my ($ctg_num) = $name =~ /Contig(\d+)\.\d+/;
	($ctg_num) = $name =~ /Contig(\d+)/ if(!defined $ctg_num);
	
    return $ctg_num;
}
my @lines = <$fh>;
chomp @lines;
my %list;
foreach (@lines)
{
	my @tokens = split;
	$list{$tokens[0]} = [$tokens[0], $tokens[1]];
}
foreach (keys %list)
{
	next if(!exists $list{$_});
	my @temp = @{$list{$_}};
	my $last_element = $temp[@temp-1];
	while(exists $list{$last_element})
	{
		shift @{$list{$last_element}};
		push @{$list{$_}}, @{$list{$last_element}};
		
		@temp = @{$list{$last_element}};
		delete $list{$last_element};
		$last_element = $temp[@temp-1];
		
	} 
}
#foreach (keys %list)
#{
#	print join " ",@{$list{$_}},"\n";
#}
#exit;
my $count = 0;

foreach (values %list)
{
	my @contigs = @{$_};
	my @ace_files = 
	my $args;
	foreach (@contigs)
	{
		$args.= $prefix.get_num($_)%$mod.".$ARGV[0] $_ ";
	}
	my $cmd = "cmt.pl $args";
	print "bsub -o mergelog$count -e errorlog$count -q short -W 5 $cmd\n";	
	$count++;	
}
