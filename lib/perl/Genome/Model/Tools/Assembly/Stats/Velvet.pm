package Genome::Model::Tools::Assembly::Stats::Velvet;

use strict;
use warnings;

use Genome;
use Cwd;
use Data::Dumper;

class Genome::Model::Tools::Assembly::Stats::Velvet {
    is => ['Genome::Model::Tools::Assembly::Stats'],
    has => [
	first_tier => {
	    type => 'Integer',
	    is_optional => 1,
	    doc => "first tier value",
	},
	second_tier => {
	    type => 'Integer',
	    is_optional => 1,
	    doc => "second tier value",
	},
	assembly_directory => {
	    type => 'Text',
	    is_optional => 1,
	    doc => "path to assembly",
	},
	major_contig_length => {
	    type => 'Integer',
	    is_optional => 1,
	    doc => "Major contig length cutoff",
	},
	out_file => { #TODO - rename this output_file
	    type => 'Text',
	    is_optional => 1,
	    is_mutable => 1,
	    doc => "Stats output file name",
	},
	no_print_to_screen => {
	    is => 'Boolean',
	    is_optional => 1,
	    default_value => 0,
	    doc => 'Prevent printing of stats to screen',
	},
	msi_assembly => {
	    is => 'Boolean',
	    is_optional => 1,
	    default_value => 0,
	    doc => 'Denote msi assemblies',
	},
	report_core_gene_survey => {
	    is => 'Boolean',
	    is_optional => 1,
	    default_value => 0,
	    doc => 'Reports core gene survey results',
	},
    ],
};

sub help_brief {
    'Run stats on velvet assemblies'
}

sub help_detail {
    return <<"EOS"
gmt assemby stats velvet --assembly-directory /gscmnt/sata910/assembly/Escherichia_coli_HMPREF9530-1.0_100416.vel
gmt assemby stats velvet --assembly-directory /gscmnt/sata910/assembly/Escherichia_coli_HMPREF9530-1.0_100416.vel --out-file stats.txt --no-print-to-screen
gmt assemby stats velvet --assembly-directory /gscmnt/sata910/assembly/Escherichia_coli_HMPREF9530-1.0_100416.vel --out-file stats.txt --first-tier 1000000 --second-tier 1200000 --major-contig-length 2000
EOS
}

sub execute {
    my $self = shift;

    my $stats; #holds streams of incoming stats text

    unless ($self->validate_assembly_out_files) {
	$self->error_message("Failed to validate assembly out files");
	return;
    }
    unless ($self->validate_velvet_assembly_files) {
	$self->error_message("Failed to validate velvet assembly files");
	return;
    }

    #SIMPLE READ STATS
    my ($s_stats, $five_k_stats, $content_stats) = $self->get_simple_read_stats();
    $stats .= $s_stats;
    print $s_stats unless $self->no_print_to_screen;

    #CONTIGUITY STATS
    my $contiguity_stats = $self->get_contiguity_stats;
    $stats .= $contiguity_stats;
    print $contiguity_stats unless $self->no_print_to_screen;

    #CONSTRAINT STATS
    my $constraint_stats = $self->get_constraint_stats();
    $stats .= $constraint_stats;
    print $constraint_stats unless $self->no_print_to_screen;

    #GENOME CONTENTS
    $stats .= $content_stats;
    print $content_stats unless $self->no_print_to_screen;

    #GENE CORE SURVEY STATS - This is optional
    if ($self->report_core_gene_survey) {
	my $core_survey = $self->get_core_gene_survey_results();
	$stats .= $core_survey;
	print $core_survey unless $self->no_print_to_screen;
    }

    #READ DEPTH STATS
    my $depth_stats = $self->get_read_depth_stats_from_afg();
    $stats .= $depth_stats;
    print $depth_stats unless $self->no_print_to_screen;

    #FIVE KB CONTIG STATS
    $stats .= $five_k_stats;
    print $five_k_stats unless $self->no_print_to_screen;

    unless ($self->out_file) {
	$self->out_file($self->assembly_directory.'/edit_dir/stats.txt');
    }

    if ($self->out_file) {
	my $out_file = $self->out_file;
	#my $fh = IO::File->new(">$out_file") || die;
	unlink $self->out_file;
	my $fh = Genome::Sys->open_file_for_writing($self->out_file) ||
	    return;
	$fh->print($stats);
	$fh->close;
    }

    print "############\n##  DONE  ##\n############\n" unless $self->no_print_to_screen;

    #returning to original dir so tests can clean up
    #chdir $dir;

    return 1;
}

1;
