package Genome::Model::Tools::Assembly::DetectMerges; 

use strict;
use warnings;

use Genome;

use FindBin;
use File::Basename;
use Carp::Assert;
use Carp;
use Cwd;
use Utility;

class Genome::Model::Tools::Assembly::DetectMerges {
    is => 'Command',
    has => [
        ace_file => { 
            is => 'String', 
            shell_args_position => 1,
            is_optional => 0,
            doc => 'the file name that we are going to detect merges from, imported assembly, or prior msi to examine',
            is_input => 1,
        },
        stats_file => {
            is => 'Number',
            shell_args_position => 2,
            is_optional => 1,
            doc => 'the name of the output stats file',
        }
    ],
    
    doc => 'identify possible merges in an assembly (usable in the merge command)'
};

my %super_contigs;
my @super_contig_names;

sub execute {
    my $self = shift;

    my $ace_file = $self->ace_file;
    my $stats_file = $self->stats_file;
    unless (defined $stats_file)
    {
        $stats_file = $ace_file.'.merge_detector_stats';
    }
    
    my $ct = Genome::Model::Tools::Pcap::ContigTools->new;
    my $ao = Genome::Model::Tools::Pcap::Ace->new(input_file => $ace_file);
    my $dir = File::Basename::dirname($ace_file);
    $dir =~ s/edit_dir//;
    my $po;
    if(-e $dir."phdball_dir/phd.ball.1") {
        $po = Genome::Model::Tools::Pcap::Phd->new(input_file => $dir."/phdball_dir/phd.ball.1");
    }
    elsif(-e $dir."/phd_dir/") {
        $po = Genome::Model::Tools::Pcap::Phd->new(input_directory => $dir."/phd_dir/");
    }
    else {
        $self->error_message("Need to either have a ../phd_dir or a phdball file named ../phdball_dir/phd.ball.1");
        return;
    }    

    my @contig_names = @{$ao->get_contig_names};

    my $outfh = IO::File->new(">$stats_file");
    $self->error_message("Error opening $stats_file for writing.\n") and return unless $outfh;
    my %params = ( statsfh => $outfh, cutoffs => { hqlength => 1000, hq_percent_identity => 90 });
    sort_supercontigs(@contig_names);
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

    print $outfh "Finished detecting merges\n";

    
    
    return 1;
}



#subroutines
sub sort_supercontigs
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
    my ($ctg_num) = $name =~ /Contig(\d+)\.\d+/;
	($ctg_num) = $name =~ /Contig(\d+)/ if(!defined $ctg_num);
	
    return $ctg_num;
}

sub get_supercontig_name
{
	my ($name) = @_;
	my ($ctg_num) = $name =~ /(Contig\d+)\.\d+/;
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

1;
