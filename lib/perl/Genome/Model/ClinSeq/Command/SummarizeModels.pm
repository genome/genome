package Genome::Model::ClinSeq::Command::SummarizeModels;
use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::SummarizeModels {
    is => 'Command::V2',
    has_input => [
        models => { 
            is => 'Genome::Model::ClinSeq',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'clinseq models to sumamrize'
        },
        outdir => { 
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written', 
        },
        skip_lims_reports => {
            is => 'Number',
            is_optional => 1,
            doc => 'Use this option to skip LIMS use of illumina_info run/lane/library report tool',
        }
    ],
    doc => 'summarize clinseq model status, input models, processing profiles, and results',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq summarize-models --outdir=/tmp/ 2882726707

genome model clin-seq summarize-models --outdir=/tmp/ "name='ClinSeq - ALL1 - (Nov. 2011 PP) - v2'"

genome model clin-seq summarize-models --outdir=/tmp/ subject.common_name=HG1

genome model clin-seq summarize-models --outdir=/tmp/ "subject.common_name like 'HG%'"

EOS
}

sub help_detail {
    return <<EOS
Summarize the status and key metrics for 1 or more clinseq models.

(put more content here)
EOS
}

sub execute {
  my $self = shift;
  my @models = $self->models;
    
  for my $model (@models) {
    #Check for a complete build of the clinseq model specified
    my $clinseq_build = $model->last_complete_build;
    unless ($clinseq_build) {
      $self->status_message("\n***** " . $model->__display_name__ . " ****");
      $self->status_message("\n\nSamples and instrument data");
      $self->status_message("NO COMPLETE CLINSEQ BUILD!");
      next;
    }
    
    my $run = Genome::Model::ClinSeq::Command::SummarizeBuilds->create(outdir=>$self->outdir, skip_lims_reports=>1, builds=>[$clinseq_build]);
    $run->dump_status_messages(1);
    $run->execute;

  }

  $self->status_message("\n\n");

  return 1;
}

1;

