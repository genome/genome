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
        dry_run             => { is => 'Boolean', doc => 'skip doing actual work (tests API)', },
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

sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    my $data_directory = $build->data_directory;

    my $wgs_build           = $build->wgs_build;
    my $exome_build         = $build->exome_build;
    my $tumor_rnaseq_build  = $build->tumor_rnaseq_build;
    my $normal_rnaseq_build = $build->normal_rnaseq_build;
    
    #This input is used for testing, and when set will not actually do any work just organize params
    my $dry_run = $build->dry_run;

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

    # initial inputs are for the "Main" component which does 
    # all of the legacy/non-parellel tasks
    my @inputs = (
        build => $build,
        wgs_som_var_data_set => $wgs_build,
        exome_som_var_data_set => $exome_build,
        tumor_rna_seq_data_set => $tumor_rnaseq_build,
        normal_rna_seq_data_set => $normal_rnaseq_build,
        working_dir => $data_directory,
        common_name => $final_name,
        verbose => 1,
        clean => 0, # unfriendly to parallelization ..dir is new and will be clean
        dry_run => $dry_run,
    );

    my $patient_dir = $data_directory . "/" . $final_name;
    my @dirs = ($patient_dir);

    # summarize builds
    my $input_summary_dir = $patient_dir . "/input";
    push @dirs, $input_summary_dir;
    push @inputs, summarize_builds_outdir => $input_summary_dir;
    push @inputs, summarize_builds_log_file => $input_summary_dir . "/SummarizeBuilds.log.tsv";

    if ($build->id > 0) {
        push @inputs, summarize_builds_skip_lims_reports => 0;
    }
    else {
        #Watch out for -ve build IDs which will occur when the ClinSeq.t test is run.  In that case, do not run the LIMS reports
        push @inputs, summarize_builds_skip_lims_reports => 1;
    }

    # For now it works to create directories here because the data_directory
    # has been allocated.  It is possible that this would not happen until
    # later, which would mess up assigning inputs to many of the commands.
    for my $dir (@dirs) {
        Genome::Sys->create_directory($input_summary_dir);
    }

    return @inputs;
}

sub _resolve_workflow_for_build {
    # This is called by Genome::Model::Build::start()
    # Returns a Workflow::Operation
    my $self = shift;
    my $build = shift;
    my $lsf_queue = shift; # TODO: the workflow shouldn't need this yet
    my $lsf_project = shift;

    if (!defined $lsf_queue || $lsf_queue eq '' || $lsf_queue eq 'inline') {
        $lsf_queue = 'apipe';
    }
    if (!defined $lsf_project || $lsf_project eq '') {
        $lsf_project = 'build' . $build->id;
    }

    # The wf system will call this method after we finish
    # but it consolidates logic to make a workflow which will
    # automatically take all inputs that method will assign.
    my %inputs = $self->map_workflow_inputs($build);
    my @input_properties = sort keys %inputs;
   
    # This must be updated for each new tool added which is "terminal" in the workflow!
    # (too bad it can't just be inferred from a dynamically expanding output connector)
    my @output_properties = qw(
        summarize_builds_result
        summarize_svs_result
        summarize_cnvs_result
    );

    my $workflow = Workflow::Model->create(
        name => $build->workflow_name,
        input_properties => \@input_properties, 
        output_properties => \@output_properties,
    );

    my $log_directory = $build->log_directory;
    $workflow->log_dir($log_directory);

    my $input_connector = $workflow->get_input_connector;
    my $output_connector = $workflow->get_output_connector;

    my %steps_by_name; 
    my $step = 0;
    my $add_step = sub {
        my ($name, $cmd) = @_;

        if (substr($cmd,0,2) eq '::') {
            $cmd = 'Genome::Model::ClinSeq::Command' . $cmd;
        }
        
        unless ($cmd->can("execute")) {
            die "bad command $cmd!";
        }

        die "$name already used!" if $steps_by_name{$name};

        my $op = $workflow->add_operation(
            name => $name,
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => $cmd,
            )
        );
        $op->operation_type->lsf_queue($lsf_queue);
        $op->operation_type->lsf_project($lsf_project);
        $steps_by_name{$name} = $op;
        return $op;
    };

    my $add_link = sub {
        my ($from_op, $from_p, $to_op, $to_p) = @_;
        my $link = $workflow->add_link(
            left_operation => $from_op,
            left_property => $from_p,
            right_operation => $to_op,
            right_property => $to_p
        );
        $link or die "Failed to make link from $link from $from_p to $to_p!";
        return $link;
    };
  


    #
    # Add steps which go in parallel with the Main step before setting up Main.
    # This will ensure testing goes more quickly because these will happen first.
    #
    
    unless ($build->dry_run) {
        #Summarize build inputs using SummarizeBuilds.pm
        $step++; 
        my $msg = "Step $step. Creating a summary of input builds using summarize-builds";
        my $summarize_builds_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeBuilds");
        $add_link->($input_connector, 'build', $summarize_builds_op, 'builds');
        $add_link->($input_connector, 'summarize_builds_outdir', $summarize_builds_op, 'outdir');
        $add_link->($input_connector, 'summarize_builds_skip_lims_reports', $summarize_builds_op, 'skip_lims_reports');
        $add_link->($input_connector, 'summarize_builds_log_file', $summarize_builds_op, 'log_file');
        $add_link->($summarize_builds_op, 'result', $output_connector, 'summarize_builds_result');
    }

    #
    # Main contains all of the original clinseq non-parallel code.
    # Re-factor out pieces which are completely independent, or which are at the end first.
    #

    $step++;
    my $main_op = $add_step->("Step $step. ClinSeq Main", "Genome::Model::ClinSeq::Command::Main");
    for my $in (qw/ 
        build
        wgs_som_var_data_set
        exome_som_var_data_set
        tumor_rna_seq_data_set
        normal_rna_seq_data_set
        working_dir
        common_name
        verbose
        clean
        dry_run
    /) {
        $add_link->($input_connector, $in, $main_op, $in);
    }
  

    #
    # When the dry_run input is set, skip the rest of the workflow (for testing).
    #

    if ($build->dry_run) {
        # skip the rest of the workflow and have the result from Main 
        # result stand-in for all of the results from other steps.
        for my $out (@output_properties) {
            # If we enable a step for use during dry-run skip
            # and it links to the final output connector, skip it here
            ## next if $out =~ /^summarize_builds/;
            $add_link->($main_op, 'result', $output_connector, $out);
        }
        # bail and return the one step workflow
        # with dry-run set it won't do any work
        return $workflow;
    }
   
    #
    # Add steps which follow the main step...
    # As we refactor start with things which do not feed into anything else down-stream,
    # or which only feed into things which have already been broken-out.
    #

    #Generate a summary of SV results from the WGS SV results
    if ($build->wgs_build) {
        $step++; 
        my $msg = "Step $step. Summarizing SV results from WGS somatic variation";
        my $summarize_svs_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeSvs");
        $add_link->($main_op, 'wgs_som_var_data_set', $summarize_svs_op, 'builds');
        $add_link->($main_op, 'sv_summary_dir', $summarize_svs_op, 'outdir');
        $add_link->($summarize_svs_op, 'result', $output_connector, 'summarize_svs_result');
    }

    #Generate a summary of CNV results, copy cnvs.hq, cnvs.png, single-bam copy number plot PDF, etc. to the cnv directory
    if ($build->wgs_build) {
        $step++;
        my $msg = "Step $step. Summarizing CNV results from WGS somatic variation";
        my $summarize_cnvs_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeCnvs");
        $add_link->($main_op, 'build', $summarize_cnvs_op, 'builds');
        $add_link->($main_op, 'cnv_summary_dir', $summarize_cnvs_op, 'outdir');
        $add_link->($summarize_cnvs_op, 'result', $output_connector, 'summarize_cnvs_result');
    }

    # REMINDER:
    # For new steps be sure to add their result to the output connector if they do not feed into another step.
    # When you do that, expand the list of output properties at line 182 above. 

    my @errors = $workflow->validate();
    if (@errors) {
        for my $error (@errors) {
            $self->error_message($error);
        }
        die "Invalid workflow!";
    }

    return $workflow;
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

