package Genome::Model::ClinSeq::Command::MakeCircosPlot;
use strict;
use warnings;
use Genome; 
use Data::Dumper;

# Written by Ben Ainscough and Scott Smith, based on prototype from Obi Griffith
# See JIRA issue https://jira.gsc.wustl.edu/browse/TD-691

class Genome::Model::ClinSeq::Command::MakeCircosPlot {
    is => 'Command::V2',
    has => [
        build                   => { is => 'Genome::Model::Build::ClinSeq',
                                     doc => 'Clinseq build' },

        output_directory        => { is => 'FilesystemPath',
                                     doc => 'Directory where output will be written', },

    ],
    doc => 'This script attempts to get read counts, frequencies and gene expression values for a series of genome positions',
};

sub sub_command_category { 'pipeline' }

sub help_detail {
  return <<EOS
Generate a Circos plot for a clinseq build.
EOS
}

sub help_synopsis {
  return <<EOS
    genome model clin-seq make-circos-plot --build model.name~bainsc%clinseq1,is_last_complete=1 --output-directory /tmp/outdir
EOS
}


sub execute {
    my $self = shift;
    $self->warning_message("!!!!!!!!!!!!!!! This module is under development and does not run yet !!!!!!!!!!!!!!!!");
   
    # grap params from $self
    my $build = $self->build;
    my $output_directory = $self->output_directory;
    $self->status_message("Running on build " . $build->__display_name__); 
    $self->status_message("Output directory is " . $self->output_directory);

    # initialize directories
    unless (-d $output_directory) {
        # this module has wrappers which do logging, throw exceptions, around regular tasks
        Genome::Sys->create_directory($output_directory);
    }

    # example accessing data:
    
    # get an input from the build
    my $wgs_build = $build->wgs_build;
    
    # all inputs are optional, so test for it being there before using it
    if ($wgs_build) {
        $self->status_message("This clinseq build has WGS somatic data, build: " . $wgs_build->__display_name__);
       
        # when you have a build, you can get the path for a file in the data directory through an API call
        my $wgs_somatic_snvs_tier1_path = $wgs_build->data_set_path('effects/snvs.hq.novel.tier1','2','bed');
        $self->status_message("WGS somatic variants are at $wgs_somatic_snvs_tier1_path");
        
        # the above is the same as the following, but the following circumvents the API
        #my $wgs_somatic_snvs_tier1_path = $wgs_build->data_directory . '/effects/snvs.hq.novel.tier1.v2.bed'; 
    }
    else {
        $self->status_message("No WGS data on this clinseq model.");
    }

    # get an input from an input
    my $wgs_tumor_refalign = $wgs_build->tumor_build;

    # eventually you do this 
    # Genome::Sys->shellcmd(cmd => "cd $output_direcotry; circos");

    return 1;
}

1;

