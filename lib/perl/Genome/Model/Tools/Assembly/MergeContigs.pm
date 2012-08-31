package Genome::Model::Tools::Assembly::MergeContigs;

#usage cmt.pl infile.ace ContigName ContigName2

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Assembly::MergeContigs
{
    is => 'Command',
    has => 
    [
        contigs => {
            type => "String",
            optional => 0,
            is_input => 1,
            doc => "This is the list of input Contigs, in the format:\n --contigs 'Contig1 Contig2 Contig6'\n it is required that an ace file is listed before each contig, even in the case where all contigs are in the same ace file"
        }, 
	    gext => {
            type => "String",
            optional => 1,
		    doc => "gap extension penalty",
            default_value => 1,
	    },
	    gopen => {
            type => "String",
            optional => 1,
		    doc => "gap open penalty",
            default_value => 1,
	    },
	    ggext => {
            type => "String",
            optional => 1,
		    doc => "global gap extension penalty",
            default_value => 15,
	    },
	    ggopen => {
            type => "String",
            optional => 1,
		    doc => "global gap open penalty",
            default_value => 15,
	    },
	    ug => {
            type => "String",
            optional => 1,
		    doc => "use global alignment",
            default_value => 0,
	    },
	    q => {
            type => "Boolean",
            optional => 1,
		    doc => "quiet mode (no stats info)",
            default_value => 0,
	    },
	    'hq_percent_identity' => {
            type => "String",
            optional => 1,
		    doc => "high quality percent identity cutoff",
	    },
	    'hq_mismatch' => {
            type => "String",
            optional => 1,
		    doc => "high quality mismatch cutoff",
	    },
	    'percent_identity' => {
            type => "String",
            optional => 1,
		    doc => "percent identity cutoff",
	    },
	    mismatch => {
            type => "String",
            optional => 1,
		    doc => "mismatch cutoff",
	    },
	    'length' => {
            type => "String",
            optional => 1,
		    doc => "alignment region length cutoff",
	    },
        input_ace_file => {
            type => "String",
            optional => 1,
		    doc => "use this param to specify a single input ace file",
	    },
        cache_dir => {
            type => "String",
            optional => 1,
            is_input => 1,
		    doc => "the cache dir used.  This is a requirement for doing merges with more than one input ace file",
	    },
        output_ace_file => {
            type => "String",
            optional => 1,
		    doc => "the single ace file to write, this is only recommened for either small assemblies, or where a single ace file was used as input",
	    },
        output_ace_dir => 
        {
            type => "String",
            optional => 1,
		    doc => "output directory to write ace files to",
	    },
        output_ace_prefix =>
        {
            type => "String",
            optional => 1,
		    doc => "prefix of each ace file in the output dir, otherwise, ace files will be named 0.ace, 1.ace, etc",        
        },
        output_file_number =>
        {
            type => "Number",
            optional => 1,
            doc => "number of output files to create"
        }
    ]
};

sub help_brief {
    ""
}

sub help_synopsis { 
    return;
}
sub help_detail {
    return <<EOS 
    merge-contigs --contigs 'contig1 contig2 contig3' 
    merge-contigs --contigs 'contig1 contig2 contig3'    
EOS
}

sub execute
{
    my $self = shift;
    my %cutoffs;
    $cutoffs{hq_percent_identity} = $self->hq_percent_identity if defined $self->hq_percent_identity;
    $cutoffs{hq_mismatch} = $self->hq_mismatch if defined $self->hq_mismatch;
    $cutoffs{percent_identity} = $self->percent_identity if defined $self->percent_identity;
    $cutoffs{mismatch} = $self->mismatch if defined $self->mismatch; 
    $cutoffs{length} = $self->length if defined $self->length;
    
    my %params = (gap_ext_penalty => $self->gext,
			      gap_open_penalty => $self->gopen,
			      glob_gap_ext_penalty => $self->ggext,
			      glob_gap_open_penalty => $self->ggopen,
			      use_global_align => $self->ug,
			      quiet => $self->q,
			      cutoffs => \%cutoffs);
    my $output_ace_dir = $self->output_ace_dir;
    my $output_ace_file = $self->output_ace_file;
    my $input_ace_file = $self->input_ace_file;
    my $cache_dir = $self->cache_dir;
    my $output_file_number = $self->output_file_number;
    my $output_ace_prefix = $self->output_ace_prefix;
    
    my $contigs = $self->contigs;
    my @contigs = split /\s+/,$contigs;
    my $ao;
    if(defined $input_ace_file)
    {
        $self->error_message("$input_ace_file is not a valid ace file, contig string $contigs is not formatted properly, there needs to be a valid ace file specified before each contig. i.e.\n merge-contigs -contigs 'contig1 contig2\n")
        and return unless (-e $input_ace_file);
        $ao = Genome::Model::Tools::Pcap::Ace->new(input_file => $input_ace_file, cache_dir => $cache_dir);
    }
    elsif(defined $cache_dir)
    {
        $self->error_message("$cache_dir is not a valid ace file, contig string $contigs is not formatted properly, there needs to be a valid ace file specified before each contig. i.e.\n merge-contigs -contigs 'contig1 contig2 '\n")
        and return unless (-e $cache_dir);
        $ao = Genome::Model::Tools::Pcap::Ace->new(cache_dir => $cache_dir);    
    }

    

    my $phd_object;
    if(-e "../phdball_dir/phd.ball.1")
    {
        $phd_object = Genome::Model::Tools::Pcap::Phd->new(input_file => "../phdball_dir/phd.ball.1");
    }
    elsif(-e "../phd_dir/")
    {
        $phd_object = Genome::Model::Tools::Pcap::Phd->new(input_directory => "../phd_dir/");
    }
    else
    {
        $self->error_message("Need to either have a ../phd_dir or a phdball file named ../phdball_dir/phd.ball.1");
        return;
    }    
    my $ct = Genome::Model::Tools::Pcap::ContigTools->new;

    my $merge_contig = $ao->get_contig($contigs[0],1);

    for(my $i=1;$i<@contigs;$i++)
    {
        my $next_contig = $ao->get_contig($contigs[$i],1);
        $merge_contig = $ct->merge($merge_contig, $next_contig, $phd_object, %params);
        $ao->remove_contig($contigs[$i]);        
    }

    $ao->add_contig($merge_contig);
    if(defined $output_ace_file)
    {
        print "Writing to output file: $output_ace_file\n";
        $ao->write_file(output_file => "$output_ace_dir/$output_ace_file");	        
    }
    elsif(defined $output_ace_dir)
    {
        `mkdir -p $output_ace_dir` if -d $output_ace_dir;
        print "Writing to output file: $output_ace_dir\n";
        $ao->write_directory(output_directory => "$output_ace_dir", prefix => $output_ace_prefix, number => $output_file_number);
    }
    
    return 1;
}
=head1 NAME

cmt : Contig Merging Tool

=head1 SYNOPSIS

cmt inputacefile leftcontig rightcontig [options]

=head1 DESCRIPTION

The cmt takes a single ace file, the name of the "left" contig and "right" contig, and produces an output ace file that
is specified by the user.  

In order for the cmt to work properly it needs quality information for the reads contained in the contigs, either
in the form of phd files.  If using phd files, they need to be in a directory named phd_dir that is one level below 
the current directory (just like consed).

	project_directory/edit_dir/(this is where are ace file are located)
	project_directory/phd_dir/(this directory contains phd information)

The cmt should be run from the directory that contains ace files, in this case edit_dir.  The directory containing ace files does not necesarily need to be called edit_dir.  However,
it is very important that the directory containing phd files is named phd_dir and the directory containing qual files is named
Input.  

=head1 OPTIONS

=head2 -gext gap-ext-penalty

Set the gap extension penalty for the local alignment algorithm.  The cmt uses a local alignment algorithm to line up the merges.  Increasing the gap extension penalty discourages the cmt from adding gaps when performing alignments.  The default value for the gap extension penalty is 1.

=head2 -gopen gap-open-penalty

Set the gap open penalty for the local alignment algorithm.   The cmt uses a local alignment algorithm to line up the merges.  Increasing the gap open penalty discourages the cmt from adding additional gaps (pads) after inserting a first gap.  The default value for gap open penalty is 1.

=head2 -ggext global-gap-open-penalty

Set the gap extension penalty for the global alignment algorithm.  The cmt uses a global alignment algorithm to line up the merges.  Increasing the gap extension penalty discourages the cmt from adding gaps when performing alignments.  The default value for global gap extension penalty is 15.

=head2 -ggopen global-gap-ext-penalty

Set the gap open penalty for the global alignment algorithm.   The cmt uses a global alignment algorithm to line parts of the Contigs that are excluded from the initial local alignment.  Increasing the gap open penalty discourages the cmt from adding additional gaps (pads) after inserting a first gap.  The default value for global gap open penalty is 15.

=head2 -ug use-global-align

Specifying this option will tell the merging tool to use global alignment to line up overlapping consensus that occurs outside the merge region.  By default this is off, since alignment of sequence outside the merge region has not been tested extensively.

=head2 -o outfilename

This option allows the user to specify the outfile name.  Otherwise the cmt will produce an output ace file that contains the input ace file plus a number.  
i.e. If our input ace file is named testin.ace then our output ace file would be testin.ace.1.

=head1 AUTHOR

The cmt was written and is maintained by Jon Schindler <jschindl@wugsc.wustl.edu>.  It is part of the contig merging and splitting toolkit.  

=head1 SEE ALSO

=cut

1;





