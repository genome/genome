package Genome::Model::Tools::Assembly::CreateMergeList; 

use IO::File;

class Genome::Model::Tools::Assembly::CreateMergeList {
    is => 'Command',
    has => [
        stats_files => {
            is => 'String',
            shell_args_position => 1,
            is_optional => 0,
            doc => 'one or more input stats files',
        },
        merge_list => {
            is => 'String',
            shell_args_position => 2,
            is_optional => 0,
            doc => 'the list of output joins',
        },
        percent_identity => {
            is => 'Number',
            default_value => 90,
            shell_args_position => 3,
            is_optional => 1,
            doc => 'the percentage of bases in the alignment region that are a match',
        },
        match_length => {
            is => 'Number',
            default_value => 200,
            shell_args_position => 4,
            is_optional => 1,
            doc => 'the length of the aligned region',
        },
    ],
    
    doc => 'Takes a  list of merge stats files produced by gmt tools assembly detect-merges and produces a list of valid joins (usable by gmt tools assembly merge-contigs)'
};



sub execute {
    my $self = shift;
    
    my $stats_files = $self->stats_files;
    my @stats_files = split (/,/,$stats_files);
    my $merge_list = $self->merge_list;
    my $percent_identity = $self->percent_identity||90;
    my $match_length = $self->match_length||200;
    my $ml_fh = IO::File->new(">$merge_list");
    $self->error_message("Could not open $merge_list for writing.\n") and return unless (defined $ml_fh);
    foreach my $file (@stats_files)
    {
	    my $fh = IO::File->new($file);
	    my $hqb = 0;
	    my $hqi = 0;
	    my $left_contig;
	    my $right_contig;
	    while(<$fh>)
	    {
		    if(/Number of high quality bases/)
		    {
			    my @tokens = split;
			    $hqb = $tokens[5];
		    }
		    if(/Merging /)
		    {
			    my @tokens = split;
			    if($hqi >= $percent_identity && $hqb >= $match_length && $left_contig && $right_contig)
			    {
				    print $ml_fh "$left_contig $right_contig $hqi $hqb\n";
			    }
			    $left_contig = $tokens[1];
			    $right_contig = $tokens[3];		
		    }
		    if(/High Quality Base Percent/)
		    {
			    my @tokens = split;
			    $hqi = $tokens[5];
		    }		
	    }
	    if($hqi >= $percent_identity && $hqb >= $match_length && $left_contig && $right_contig)
	    {
		    print $ml_fh "$left_contig $right_contig $hqi $hqb\n";
	    }	
    }
    return 1;
}

1;
