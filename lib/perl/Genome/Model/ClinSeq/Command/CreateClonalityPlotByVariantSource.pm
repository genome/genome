package Genome::Model::ClinSeq::Command::CreateClonalityPlotByVariantSource;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::CreateClonalityPlotByVariantSource {
    is => 'Command::V2',
    has_input => [
        builds => { 
              is => 'Genome::Model::Build::ClinSeq',
              is_many => 1,
              shell_args_position => 1,
              require_user_verify => 0,
              doc => 'somatic variation build(s) to get variant sources from',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
        callers => {
            is => 'String',
            doc => 'Names of the callers to graph',
            is_many => 1,
            valid_values => [qw(sniper strelka varscan)],
        },
    ],
    has_output => [
        _plots => {
              is => 'FilesystemPath',
              is_many => 1,
              is_optional =>1,
        },
    ],
    doc => 'plot clonality graphs for each of the sources of variants (i.e., which variant callers) for a ClinSeq build',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq create-clonality-plot-by-variant-source --outdir=/tmp/  --callers strelka,sniper 128884819

EOS
}

sub help_detail {
    return <<EOS
Graph Clonality plots for each combination of source of variants (i.e., snv/indel caller) for ClinSeq build

(put more content here)
EOS
}

sub execute {
    my $self = shift;
    #For each clinseq build
    for my $build ($self->builds) {
    #grab the variant source file
    #the readcount files
    #the CNAseq files
    #
    #Then, generate per-caller readcount files and their unions

    #TODO This should really operate based on DV2 results and combinations, but
    #this is more informative in the short term and the variant sources code
    #already suffers from the same issues 
    #
    #Then plot various graphs
    }
}

sub _clinseq_patient_dir {
    my ($self, $clinseq_build) = @_;
    my $patient_dir = $clinseq_build->data_directory . "/" . $clinseq_build->common_name;
    unless(-d $patient_dir) {
        $self->error_message("ClinSeq patient directory not found. Expected: $patient_dir");
        die;
    }
    return $patient_dir;
}

sub _snv_variant_source_file {
    my ($self, $clinseq_build, $data_type,) = @_;

    my $patient_dir = $self->_clinseq_patient_dir($clinseq_build);
    my $source = $patient_dir . "/variant_source_callers/";
    my $dir;
    if(-d $source) {
       $dir = $source . "/$data_type/";
       if(-d $dir) {
           my $file = $dir . "/snv_sources.tsv";
           if(-e $file) {
               return $file;
           }
           else {
               $self->error_message("Expected $file inside $dir and it did not exist.");
           }
       }
       else {
           $self->error_message("$data_type sub-directory not found in $source."); 
       }
    }
    else {
        $self->error_message("$source directory not found");
    }
    die;
}

sub _clinseq_clonality_dir {
    my ($self, $clinseq_build) = @_;

    my $patient_dir = $self->_clinseq_patient_dir($clinseq_build);
    my $clonality_dir = $patient_dir . "/clonality/";
    
    unless(-d $clonality_dir) {
        $self->error_message("Clonality directory does not exist. Expected: $clonality_dir");
        die;
    }
    return $clonality_dir;
}

sub _varscan_formatted_readcount_file {
    my ($self, $clinseq_build) = @_;
    my $clonality_dir = $self->clonality_dir($clinseq_build);
    my $readcount_file = $clonality_dir ."/allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts.varscan";
    unless(-e $readcount_file) {
        $self->error_message("Unable to find varscan formatted readcount file. Expected: $readcount_file");
        die;
    }
    return $readcount_file;
}

sub _cnaseq_hmm_file {
    my ($self, $clinseq_build) = @_;
    my $clonality_dir = $self->clonality_dir($clinseq_build);
    my $hmm_file = $clonality_dir . "/cnaseq.cnvhmm";
    unless(-e $hmm_file) {
        $self->error_message("Unable to find cnaseq hmm file. Expected: $hmm_file");
        die;
    }
    return $hmm_file;
}

sub _create_readcount_file_for_caller {
    my ($self, $reacount_file, $caller) = @_;
    
}



1;


