package Genome::Model::ClinSeq;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq {
    is => 'Genome::ModelDeprecated',
    has_optional_input => [
        wgs_model           => { is => 'Genome::Model::SomaticVariation', doc => 'somatic variation model for wgs data' },
        exome_model         => { is => 'Genome::Model::SomaticVariation', doc => 'somatic variation model for exome data' },
        tumor_rnaseq_model  => { is => 'Genome::Model::RnaSeq', doc => 'rnaseq model for tumor rna-seq data' },
        normal_rnaseq_model => { is => 'Genome::Model::RnaSeq', doc => 'rnaseq model for normal rna-seq data' },
        force               => { is => 'Boolean', doc => 'skip sanity checks on input models' },  
    ],
    has_optional_param => [
        #Processing profile parameters would go in here
        #someparam1 => { is => 'Number', doc => 'blah' },
        #someparam2 => { is => 'Boolean', doc => 'blah' },
        #someparam2 => { is => 'Text', valid_values => ['a','b','c'], doc => 'blah' },
    ],
    doc => 'clinial sequencing data convergence of RNASeq, WGS and exome capture data',
};

sub _help_synopsis {
    my $self = shift;
    return <<"EOS"

    genome processing-profile create clin-seq  --name 'November 2011 Clinical Sequencing' 

    genome model define clin-seq  --processing-profile='November 2011 Clinical Sequencing'  --wgs-model='2882504846'  --exome-model='2882505032'  --tumor-rnaseq-model='2880794613'
    
    # Automatically builds if/when the models have a complete underlying build
EOS
}

sub _help_detail_for_profile_create {
    return <<EOS

The initial ClinSeq pipeline has no parameters.  Just use the default profile to run it.

EOS
}

sub _help_detail_for_model_define {
    return <<EOS

The ClinSeq pipeline takes four models, each of which is optional, and produces data sets potentially useful in a clinical setting.

EOS
}

#TODO:  Is the following test actually performed?  I don't see it being called anywhere?
sub _resolve_subject {
    my $self = shift;
    my @subjects = $self->_infer_candidate_subjects_from_input_models();
    if (@subjects > 1) {
      if ($self->force){
        @subjects = ($subjects[0]);
      }else{
        $self->error_message(
            "Conflicting subjects on input models!:\n\t"
            . join("\n\t", map { $_->__display_name__ } @subjects)
        );
        return;
      }
    }
    elsif (@subjects == 0) {
        $self->error_message("No subjects on input models?  Contact Informatics.");
        return;
    }
    return $subjects[0];
}

#TODO: As above, infer reference genome builds from input models and abort if they are missing or in conflict!

#TODO: As above, infer annotation builds from input models and abort if they are missing or in conflict!

#Implement specific error checking here, any error that is added to the @errors array will prevent the model from being commited to the database
#Could also implement a --force input above to allow over-riding of errors
sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;

#    unless(-e $self->input_fasta_file && -f $self->input_fasta_file) {
#        push @errors,UR::Object::Tag->create(
#            type => 'error',
#            properties => ['input_fasta_file'],
#            desc => 'input_fasta_file does not exist or is not a file'
#        );
#    }

    return @errors;
}


sub _execute_build {
    my ($self,$build) = @_;

    #Make sure the right libraries are used (in case someone runs with a perl -I statement).
    my $prefix = UR::Util->used_libs_perl5lib_prefix;
    local $ENV{PERL5LIB} = $prefix . ':' . $ENV{PERL5LIB};

    my $data_directory = $build->data_directory;

    my $wgs_build           = $build->wgs_build;
    my $exome_build         = $build->exome_build;
    my $tumor_rnaseq_build  = $build->tumor_rnaseq_build;
    my $normal_rnaseq_build = $build->normal_rnaseq_build;
    
    #This input is used for testing, and when set will not actually do any work just organize params
    my $dry_run = $build->inputs(name => 'dry_run');
    $dry_run = $dry_run->value if $dry_run;

    #Note that the option names are being truncated below here deliberately.  Getopt will allow the names to be arbitrarily shorter as long as they are still unique/unambiguous
    my ($wgs_common_name, $exome_common_name, $tumor_rnaseq_common_name, $normal_rnaseq_common_name, $wgs_name, $exome_name, $tumor_rnaseq_name, $normal_rnaseq_name);
    if ($wgs_build) {
        $wgs_common_name = $wgs_build->subject->patient->common_name;
        $wgs_name = $wgs_build->subject->patient->name;
    }
    if ($exome_build) {
        $exome_common_name = $exome_build->subject->patient->common_name;
        $exome_name = $exome_build->subject->patient->name;
    }
    if ($tumor_rnaseq_build) {
        $tumor_rnaseq_common_name = $tumor_rnaseq_build->subject->patient->common_name;
        $tumor_rnaseq_name = $tumor_rnaseq_build->subject->patient->name;
    }
    if ($normal_rnaseq_build) {
        $normal_rnaseq_common_name = $normal_rnaseq_build->subject->patient->common_name;
        $normal_rnaseq_name = $normal_rnaseq_build->subject->patient->name;
    }

    #Get the patient common name from one of the builds, if none can be found, use the individual name instead, if that can't be found either set the name to 'UnknownName'
    my @names = ($wgs_common_name, $exome_common_name, $tumor_rnaseq_common_name, $normal_rnaseq_common_name, $wgs_name, $exome_name, $tumor_rnaseq_name, $normal_rnaseq_name);
    my $final_name = "UnknownName";
    foreach my $name (@names){
      if ($name){
        $final_name = $name;
        last();
      }
    }

    #Before executing, change the environment variable for R_LIBS to be ''
    #This will force R to use it own local notion of library paths instead of the /gsc/ versions
    #This should work for R installed on the machine /usr/bin/R  OR  a standalone version of R installed by a local user. e.g. /gscmnt/gc2142/techd/tools/R/R-2.14.0/bin/R
    local $ENV{R_LIBS}='';

    my $cmd = Genome::Model::ClinSeq::Command::Main->create(
      build_id => $build->id,
      ($wgs_build ? (wgs_som_var_data_set => $wgs_build->id) : ()),
      ($exome_build ? (exome_som_var_data_set => $exome_build->id) : ()),
      ($tumor_rnaseq_build ? (tumor_rna_seq_data_set => $tumor_rnaseq_build->id) : ()),
      ($normal_rnaseq_build ? (normal_rna_seq_data_set => $normal_rnaseq_build->id) : ()),
      working_dir => $data_directory,
      common_name => $final_name,
      verbose => 1,
      clean => 1,
    );

    #Unless a dry run, execute the ClinSeq Main command
    if ($dry_run) {
        $build->status_message("NOT running!");
    }else {
      open(OLD, ">&STDOUT"); #Save stdout
      open(STDOUT,">&STDERR"); #Redirect stdout to go to stderr
      my $result = $cmd->execute();
      open(STDOUT, ">&OLD"); #Return stdout to its usual state
      if (!$result){
        die "Bad return value from ClinSeq::Command::Main";
      } 
      close(OLD);
    }

    return 1;
}

sub _infer_candidate_subjects_from_input_models {
    my $self = shift;
    my %subjects;
    for my $input_model (
        $self->wgs_model,
        $self->exome_model,
        $self->tumor_rnaseq_model,
        $self->normal_rnaseq_model,
    ) {
        next unless $input_model;
        my $patient;
        if ($input_model->subject->isa("Genome::Individual")) {
            $patient = $input_model->subject;
        }
        else {
            $patient = $input_model->subject->patient;
        }
        $subjects{ $patient->id } = $patient;

        # this will only work when the subject is an original tissue
        next;

        my $tumor_model;
        if ($input_model->can("tumor_model")) {
            $tumor_model = $input_model->tumor_model;
        }
        else {
            $tumor_model = $input_model;
        }
        $subjects{ $tumor_model->subject_id } = $tumor_model->subject;
    }
    my @subjects = sort { $a->id cmp $b->id } values %subjects;
    return @subjects;
}

sub _resolve_resource_requirements_for_build {
  #Set LSF parameters for the ClinSeq job
  my $lsf_resource_string = "-R 'select[model!=Opteron250 && type==LINUX64] rusage[tmp=10000:mem=4000]' -M 4000000";
  return($lsf_resource_string);
}

# This is implemented here until refactoring is done on the model/build API
# It ensures that when the CI server compares two clinseq builds it ignores the correct files.
# Keep it in sync with the diff conditions in ClinSeq.t.
sub files_ignored_by_build_diff {
    return qw(
        build.xml
        reports/Build_Initialized/report.xml
        reports/Build_Succeeded/report.xml
        logs/.*
        .*.R$
        .*.pdf$
        .*._COSMIC.svg$
        .*.clustered.data.tsv$
        .*.SummarizeBuilds.log.tsv$
        .*.DumpIgvXml.log.txt
        .*/mutation_diagrams/cosmic.mutation-diagram.stderr
        .*/mutation_diagrams/somatic.mutation-diagram.stderr
    );
};

1;

