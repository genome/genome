package Genome::Model::Tools::Analysis::PlotMultiMutationSpectrum;

#####################################################################################################################################
# MutationSpectrum - Given an annotation file, gives an output of transition/transversion, cpg islands, and cpgs within cpg islands.
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	3/05/2010 by W.S.
#	MODIFIED:	3/05/2010 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use warnings;
use strict;

use Genome;
use Workflow;
use Carp;
use FileHandle;
use Data::Dumper;
use List::Util qw( max );
use IO::File;
use Genome::Info::IUB;
use DBI;
use Cwd ('getcwd','abs_path');

class Genome::Model::Tools::Analysis::PlotMultiMutationSpectrum {
    is => ['Command'],
    has => [
    mut_spec_files => { 
        is  => 'String',
        is_input=>1, 
	is_optional => 1,
        doc => 'comma-separated mutation spectrum outputfiles',
    },
    labels => { 
        is  => 'String',
        is_input=>1, 
	is_optional => 1,
        doc => 'comma-separated labels that correspond to the inputfiles',
    },
    grouping_file => {
	is  => 'String',
        is_input=>1, 
	is_optional => 1,
        doc => 'tab-delimited file containing mapping of mutation_spectrum files and labels, supercedes --mut-spec-files and --labels',
    },
    plot_title => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => 'Mutation Spectrums',
        doc => 'The title of the plot',
    },
    plot_type => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => 'facet1',
        doc => 'Affect the type of plots made',
    },
    calc_pvalue => {
        is_input => 1,
        is_optional => 1,
	is => 'Boolean',
        default => 0,
        doc => 'Calculate pvalue for each mutation type (default is FALSE)  Will only work for plot-type=facet1.',
    },
    plot_output_file => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 0,
        #default_value => 'output.pdf',
        doc => 'The name of pdf file to save the plot to',
    },
    plot_order => {
        is_input => 1,
        is_optional => 1,
	is => 'String',
        default => 0,
        doc => 'order at which the labels will be plotted on the graph (comma,separated), will use the order defined by --label as default',
    },
    temp_file => {
	is_input => 1,
        is_optional => 1,
	is => 'String',
	doc => 'The name of the temporary file to save the data to be plotted',
    #below here are variables with which to store results
    #hopefully these can then be judiciously used to write cross-comparison scripts
    },
    _total_lines => {
        is => 'Integer',
        is_optional => 1,
        doc => "Total number of lines in the file",
    },   
    _num_indels => {
        is => 'Integer',
        is_optional => 1,
        doc => "Total number of indels in the file",
    },   
    _num_mnps => {
        is => 'Integer',
        is_optional => 1,
        doc => "Total number of mnp (likely dnps) in the file",
    },   
    _num_snvs => {
        is => 'Integer',
        is_optional => 1,
        doc => "Total number of snvs in the file",
    },   
    _transitions_transversion_count => {
        is => 'HashRef',
        is_optional => 1,
        doc => "Hash containing the counts for the various base changes as hash{base1}{base2}",
    },   
    _num_transitions => {
        is => 'Integer',
        is_optional => 1,
        doc => "Number of transitions",
    },
    _num_transversions => {
        is => 'Integer',
        is_optional => 1,
        doc => "Number of transversions",
    },

    _synonymous_transitions_transversion_count => {
        is => 'HashRef',
        is_optional => 1,
        doc => "Hash containing the counts for the various base changes for synonymous changes as hash{base1}{base2}",
    },   
    #not dealing with CpG/NpG stuff yet

    # Make workflow choose 64 bit blades
    lsf_resource => {
        is_param => 1,
        default_value => 'rusage[mem=4000] select[type==LINUX64] span[hosts=1]',
    },
    lsf_queue => {
        is_param => 1,
        default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
    }, 
    ],
};

sub help_brief {
    "Given an annotation file, gives an output of transition/transversion, cpg islands, and cpgs within cpg islands.";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"

gmt analysis plot-multi-mutation-spectrum --grouping-file luc_file  --plot-title 'Lung Cancer' --calc-pvalue --plot-output-file LUC.pdf --plot-type 'facet1'
gmt analysis plot-multi-mutation-spectrum --mut-spec-files luc1.mut.spec,luc4.mut.spec --labels luc1,luc4 --plot-title 'Lung Cancer' --calc-pvalue --plot-output-file LUC.pdf --plot-type 'facet1'
EXAMPLE:
gmt analysis plot-multi-mutation-spectrum --grouping-file luc_file  --plot-title 'Lung Cancer' --calc-pvalue --plot-output-file LUC.pdf --plot-type 'facet1'

EOS
}

sub help_detail {                           
    return <<EOS 
    Given multiple output of mutation spectrum outputfile, make plots that compares the mutation spectrum of multiple samples.  You can specify the files using either a file via --group-file (tab-delimited where 1st column is the path to the file, 2nd column is the label) or via a comma-separated list via --mut-spec-files and --labels (make sure the order and number match)..  This tool can make 4 different plots.  --plot-type=facet1 (default) separates the plots into multiplot by mutation type.  --plot-type=facet2 separates the plot into multiple plots by sample.  --plot-type=bar1 is a standard bargraph (similar to facet1 (but only makes 1 graph).  --plot-type=bar2 is a stacked barplot which is great if you have large number of samples.  
EOS
    }

sub execute {
    my $self = shift;
    $DB::single = 1;

    my $gr_file = $self->grouping_file;
    my $mutation_spec_files=[]; #arrayref
    my $input_labels=[]; #arrayref
    if(defined($gr_file)) {
	($mutation_spec_files, $input_labels) = parse_grouping_file($gr_file);
    }else {
	@$mutation_spec_files = split(",",$self->mut_spec_files);
	@$input_labels = split(",",$self->labels);
    }
    my $plot_label_order = $self->plot_order || join(",",@$input_labels);

    #my @mutation_spec_files = split(",",$self->mut_spec_files);
    #my @input_labels = split(",",$self->labels);

    if(scalar(@$mutation_spec_files) != scalar(@$input_labels)) {
	die "Number of Input Files Must Match Number of Labels!!!\n";
    }

    my $temp_file = $self->temp_file;
    my $input_plot_file;
    my ($fh,$tfile);
    if(defined($temp_file)) {
	$temp_file = abs_path($temp_file);
	unlink($temp_file) if(-e $temp_file);
	$fh = Genome::Sys->open_file_for_writing($temp_file);
	$input_plot_file = $temp_file;
    }
    else {
	($fh,$tfile) = Genome::Sys->create_temp_file;	
	$input_plot_file = $tfile;
    }



    my $header = "Category\tCount\tDensity\tSample";
    $fh->print("$header\n");

    for(my $i=0;$i<@$mutation_spec_files;$i++) {
	my $mut_spec = parse_mut_spec_file($mutation_spec_files->[$i]);
	my $lab = $input_labels->[$i];
	
	foreach my $cat(sort keys %{$mut_spec->{'basechg'}}) {
	    my $count = $mut_spec->{'basechg'}->{$cat}->[0];
	    my $density = $mut_spec->{'basechg'}->{$cat}->[1];
	    $fh->print("$cat\t$count\t$density\t$lab\n");
	}
	
	foreach my $cat(sort keys %{$mut_spec->{'class'}}) {
	    my $count = $mut_spec->{'class'}->{$cat}->[0];
	    my $density = $mut_spec->{'class'}->{$cat}->[1];
	    $fh->print("$cat\t$count\t$density\t$lab\n");
	}
    }
    close $fh;



    #my $input_plot_file = $out1;
    my $plot_output_file = $self->plot_output_file;

    if($plot_output_file !~ /^\//) { #if user does not specify absolute path, use current directory
	my $cwd = getcwd();
	$plot_output_file = "$cwd/$plot_output_file";
    }
    unless($plot_output_file =~ /\.pdf$/) {
	$plot_output_file .= ".pdf";
    }

    my $plot_title = $self->plot_title;
    my $plot_type = $self->plot_type;

    my $calc_pvalue = 'FALSE';
    if($self->calc_pvalue) {
	$calc_pvalue='TRUE'
    }

    my $plot_cmd;
    if($plot_label_order) {
	$plot_cmd = qq{ plot_multi_mutation_spectrum("$input_plot_file",output_file="$plot_output_file",plot_title="$plot_title",plot_type="$plot_type",pvalue=$calc_pvalue,plot.sample.order='$plot_label_order') };
    }

    my $call = Genome::Model::Tools::R::CallR->create(command=>$plot_cmd, library=> "MutationSpectrum.R");
    $call->execute;

    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

sub parse_mut_spec_file {

    my $file = shift;

    open(FILE, $file) or die "Can't open the file $file due to $!";
    my $data={};
    while(<FILE>) {
	chomp;
	my @list = split(/\t/,$_);
	next if($list[0] =~ /base->base_change/i);
	next if($list[0] =~ /Indels/i);
	next if($list[0] =~ /DNP0/i);

	if($list[0] =~ /->/) {
	    $data->{'basechg'}->{$list[0]} = [$list[1]];
	}
	if($list[0] =~ /Tran/i){
	    $data->{'class'}->{$list[0]} = [$list[1]];
	}
	if($list[0] =~ /SNVs/i) {
	    $data->{'total'} = $list[1];
	}
    }
    close FILE;

    foreach (keys %{$data->{'basechg'}}) {
	my $percent = ($data->{'basechg'}->{$_}->[0]/$data->{'total'})*100;
	$data->{'basechg'}->{$_}->[1] = $percent
    }
    foreach (keys %{$data->{'class'}}) {
	my $percent = ($data->{'class'}->{$_}->[0]/$data->{'total'})*100;
	$data->{'class'}->{$_}->[1] = $percent
    }

    return $data;

}


sub parse_grouping_file {

    my $file = shift;

    open(FILE, $file) or die "Can't open the file $file due to $!";
    
    my @files=();
    my @labels=();
    while(<FILE>) {
	chomp;
	next if(/^\s+$/); #remove empty lines
	my ($mut_spec_file,$label) = split(/\t/,$_);
	next if(!$mut_spec_file or !$label);
	$mut_spec_file = rem_white_space($mut_spec_file);
	$label =  rem_white_space($label);

	push(@files,$mut_spec_file);
	push(@labels,$label);
    }
    close FILE;

    return (\@files,\@labels);


}


sub rem_white_space {

    my $string = shift;

    $string =~ s/^\s+//;
    $string =~ s/\s+$//;

    return $string;
}
