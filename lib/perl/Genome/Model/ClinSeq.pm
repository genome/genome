package Genome::Model::ClinSeq;

use strict;
use warnings;
use Genome;

# these are used below, and are also used in the documentation on commands in the tree
# to provide the most useful examples possible
our $DEFAULT_CANCER_ANNOTATION_DB_ID    = 'tgi/cancer-annotation/human/build37-20130401.1';
our $DEFAULT_MISC_ANNOTATION_DB_ID      = 'tgi/misc-annotation/human/build37-20130113.1';
our $DEFAULT_COSMIC_ANNOTATION_DB_ID    = 'cosmic/65.1';

class Genome::Model::ClinSeq {
    is => 'Genome::Model',
    has_optional_input => [
        wgs_model               => { is => 'Genome::Model::SomaticVariation', doc => 'somatic variation model for wgs data' },
        exome_model             => { is => 'Genome::Model::SomaticVariation', doc => 'somatic variation model for exome data' },
        tumor_rnaseq_model      => { is => 'Genome::Model::RnaSeq', doc => 'rnaseq model for tumor rna-seq data' },
        normal_rnaseq_model     => { is => 'Genome::Model::RnaSeq', doc => 'rnaseq model for normal rna-seq data' },
        de_model                => { is => 'Genome::Model::DifferentialExpression', doc => 'differential-expression for tumor vs normal rna-seq data' },

        cancer_annotation_db    => { is => 'Genome::Db::Tgi::CancerAnnotation', default_value => $DEFAULT_CANCER_ANNOTATION_DB_ID }, 
        misc_annotation_db      => { is => 'Genome::Db::Tgi::MiscAnnotation', default_value => $DEFAULT_MISC_ANNOTATION_DB_ID },
        cosmic_annotation_db    => { is => 'Genome::Db::Cosmic', default_value => $DEFAULT_COSMIC_ANNOTATION_DB_ID },
        
        force                   => { is => 'Boolean', doc => 'skip sanity checks on input models' },

        #processing_profile      => { is => 'Genome::ProcessingProfile::ClinSeq', id_by => 'processing_profile_id', default_value => { } },
    ],
    has_optional_param => [
        #Processing profile parameters would go in here
        #someparam1 => { is => 'Number', doc => 'blah' },
        #someparam2 => { is => 'Boolean', doc => 'blah' },
        #someparam2 => { is => 'Text', valid_values => ['a','b','c'], doc => 'blah' },
    ],
    has_optional_metric => [
        common_name         => { is => 'Text', doc => 'the name chosen for the root directory in the build' },
    ],
    has_calculated => [
        expected_common_name => {
            is => 'Text',
            calculate_from => [qw /wgs_model exome_model tumor_rnaseq_model normal_rnaseq_model/],
            calculate => q|
              my ($wgs_common_name, $exome_common_name, $tumor_rnaseq_common_name, $normal_rnaseq_common_name, $wgs_name, $exome_name, $tumor_rnaseq_name, $normal_rnaseq_name);
              if ($wgs_model) {
                  $wgs_common_name = $wgs_model->subject->patient->common_name;
                  $wgs_name = $wgs_model->subject->patient->name;
              }
              if ($exome_model) {
                  $exome_common_name = $exome_model->subject->patient->common_name;
                  $exome_name = $exome_model->subject->patient->name;
              }
              if ($tumor_rnaseq_model) {
                  $tumor_rnaseq_common_name = $tumor_rnaseq_model->subject->patient->common_name;
                  $tumor_rnaseq_name = $tumor_rnaseq_model->subject->patient->name;
              }
              if ($normal_rnaseq_model) {
                  $normal_rnaseq_common_name = $normal_rnaseq_model->subject->patient->common_name;
                  $normal_rnaseq_name = $normal_rnaseq_model->subject->patient->name;
              }
              my @names = ($wgs_common_name, $exome_common_name, $tumor_rnaseq_common_name, $normal_rnaseq_common_name, $wgs_name, $exome_name, $tumor_rnaseq_name, $normal_rnaseq_name);
              my $final_name = "UnknownName";
              foreach my $name (@names){
                if ($name){
                  $final_name = $name;
                  last();
                }
              }
              return $final_name;
            |
        },
    ],
    has => [
        individual_common_name => {
            via => '__self__',
            to => 'expected_common_name',
            doc => 'alias for benefit of Solr indexing',
        },
    ],
    doc => 'clinical and discovery sequencing data analysis and convergence of RNASeq, WGS and exome capture data',
};

sub define_by { 'Genome::Model::Command::Define::BaseMinimal' }

sub _help_synopsis {
    my $self = shift;
    return <<"EOS"

    genome processing-profile create clin-seq  --name 'November 2011 Clinical Sequencing' 

    genome model define clin-seq  --processing-profile='November 2011 Clinical Sequencing'  --wgs-model=2882504846 --exome-model=2882505032 --tumor-rnaseq-model=2880794613
    
    # Automatically builds if/when the models have a complete underlying build
EOS
}

sub _help_detail_for_profile_create {
    return <<EOS

The ClinSeq pipeline has no parameters.  Just use the default profile to run it.

EOS
}

sub _help_detail_for_model_define {
    return <<EOS

The ClinSeq pipeline takes four models, each of which is optional, and produces data sets potentially useful in a clinical or discovery project.

There are several primary goals:

1. Generate results that help to identify clinically actionable events (gain of function mutations, amplifications, etc.)

2. Converge results across data types for a single case (e.g. variant read counts from wgs + exome + rnaseq)

3. Automatically generate statistics, tables and figures to be used in manuscripts

4. Test novel analysis methods in a pipeline that runs faster and is more agile than say somatic-variation

EOS
}

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
  return @errors;
}

sub _initialize_build {
    my $self = shift;
    my $build = shift;

    # this is currently tracked as a metric on the build but it should really be an output/metric
    my $common_name = $self->expected_common_name;
    $build->common_name($common_name);

    return 1;
}

#MAP WORKFLOW INPUTS
sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    my $data_directory = $build->data_directory;

    # inputs
    my $wgs_build           = $build->wgs_build;
    my $exome_build         = $build->exome_build;
    my $tumor_rnaseq_build  = $build->tumor_rnaseq_build;
    my $normal_rnaseq_build = $build->normal_rnaseq_build;
    
    # set during build initialization
    my $common_name         = $build->common_name;
    unless ($common_name) {
        die "no common name set on build during initialization?";
    }

    my $cancer_annotation_db = $build->cancer_annotation_db;
    my $misc_annotation_db = $build->misc_annotation_db;
    my $cosmic_annotation_db = $build->cosmic_annotation_db;

    # initial inputs used for various steps
    my @inputs = (
        build => $build,
        build_as_array => [$build],
        wgs_build => $wgs_build,
        wgs_build_as_array => [$wgs_build],
        exome_build => $exome_build,
        exome_build_as_array => [$exome_build],
        tumor_rnaseq_build => $tumor_rnaseq_build,
        normal_rnaseq_build => $normal_rnaseq_build,
        working_dir => $data_directory,
        common_name => $common_name,
        verbose => 1,
        cancer_annotation_db => $cancer_annotation_db,
        misc_annotation_db => $misc_annotation_db,
        cosmic_annotation_db => $cosmic_annotation_db,
    );

    my $patient_dir = $data_directory . "/" . $common_name;
    my @dirs = ($patient_dir);

    #SummarizeBuilds
    my $input_summary_dir = $patient_dir . "/input";
    push @dirs, $input_summary_dir;
    push @inputs, summarize_builds_outdir => $input_summary_dir;
    push @inputs, summarize_builds_log_file => $input_summary_dir . "/SummarizeBuilds.log.tsv";

    if ($build->id) {
      push @inputs, summarize_builds_skip_lims_reports => 0;
    }else {
      #Watch out for -ve build IDs which will occur when the ClinSeq.t test is run.  In that case, do not run the LIMS reports
      push @inputs, summarize_builds_skip_lims_reports => 1;
    }

    #DumpIgvXml
    my $igv_session_dir = $patient_dir . '/igv';
    push @dirs, $igv_session_dir;
    push @inputs, igv_session_dir => $igv_session_dir;

    #GetVariantSources
    if ($wgs_build or $exome_build) {
      my $variant_sources_dir = $patient_dir . '/variant_source_callers';
      push @dirs, $variant_sources_dir;
      if ($wgs_build){
        my $wgs_variant_sources_dir = $variant_sources_dir . '/wgs';
        push @inputs, wgs_variant_sources_dir => $wgs_variant_sources_dir;
        push @dirs, $wgs_variant_sources_dir;
      }
      if ($exome_build){
        my $exome_variant_sources_dir = $variant_sources_dir . '/exome';
        push @inputs, exome_variant_sources_dir => $exome_variant_sources_dir;
        push @dirs, $exome_variant_sources_dir;
      }
    }

    #ImportSnvsIndels
    push @inputs, import_snvs_indels_filter_mt => 1;
    push @inputs, import_snvs_indels_outdir => $patient_dir;

    #CreateMutationDiagrams
    if ($wgs_build or $exome_build) {
      my $mutation_diagram_dir = $patient_dir . '/mutation_diagrams';
      push @dirs, $mutation_diagram_dir;
      push @inputs, (
          mutation_diagram_outdir => $mutation_diagram_dir,
          mutation_diagram_collapse_variants=>1, 
          mutation_diagram_max_snvs_per_file=>1500, 
          mutation_diagram_max_indels_per_file=>1500,
          #mutation_diagram_cosmic_version=> $build->cosmic_annotation_db->external_version,
      );
    }

    #rnaseq analysis steps according to what rna-seq builds are defined as inputs
    if ($normal_rnaseq_build or $tumor_rnaseq_build){
      my $rnaseq_dir = $patient_dir . '/rnaseq';
      push @dirs, $rnaseq_dir;
      push @inputs, cufflinks_percent_cutoff => 1;

      #TophatJunctionsAbsolute and CufflinksExpressionAbsolute for 'normal' sample
      if ($normal_rnaseq_build){
        my $normal_rnaseq_dir = $rnaseq_dir . '/normal';
        push @dirs, $normal_rnaseq_dir;
        my $normal_tophat_junctions_absolute_dir = $normal_rnaseq_dir . '/tophat_junctions_absolute';
        push @dirs, $normal_tophat_junctions_absolute_dir;
        push @inputs, normal_tophat_junctions_absolute_dir => $normal_tophat_junctions_absolute_dir;
        my $normal_cufflinks_expression_absolute_dir = $normal_rnaseq_dir . '/cufflinks_expression_absolute';
        push @dirs, $normal_cufflinks_expression_absolute_dir;
        push @inputs, normal_cufflinks_expression_absolute_dir => $normal_cufflinks_expression_absolute_dir;
      }

      #TophatJunctionsAbsolute and CufflinksExpressionAbsolute for 'tumor' sample
      if ($tumor_rnaseq_build){
        my $tumor_rnaseq_dir = $rnaseq_dir . '/tumor';
        push @dirs, $tumor_rnaseq_dir;
        my $tumor_tophat_junctions_absolute_dir = $tumor_rnaseq_dir . '/tophat_junctions_absolute';
        push @dirs, $tumor_tophat_junctions_absolute_dir;
        push @inputs, tumor_tophat_junctions_absolute_dir => $tumor_tophat_junctions_absolute_dir;
        my $tumor_cufflinks_expression_absolute_dir = $tumor_rnaseq_dir . '/cufflinks_expression_absolute';
        push @dirs, $tumor_cufflinks_expression_absolute_dir;
        push @inputs, tumor_cufflinks_expression_absolute_dir => $tumor_cufflinks_expression_absolute_dir;
      }

      #CufflinksDifferentialExpression
      if ($normal_rnaseq_build and $tumor_rnaseq_build){
        my $cufflinks_differential_expression_dir = $rnaseq_dir . '/cufflinks_differential_expression';
        push @dirs, $cufflinks_differential_expression_dir;
        push @inputs, cufflinks_differential_expression_dir => $cufflinks_differential_expression_dir;
      }

      #Filtered and Intersected ChimeraScan fusion output for tumor RNAseq
      #The following will only work if the rna-seq build used a processing profile that uses chimerascan.
      if ($tumor_rnaseq_build){
        #Check for ChimeraScan fusion results
        if(-e $tumor_rnaseq_build->data_directory . '/fusions/chimeras.bedpe'){
            my $ncbi_human_ensembl_build_id = $build->tumor_rnaseq_build->annotation_build->id;
            my $tumor_unfiltered_fusion_file = $build->tumor_rnaseq_build->data_directory . '/fusions/chimeras.bedpe';
            my $tumor_filtered_fusion_dir = $patient_dir . '/rna_fusions/tumor';
            my $tumor_filtered_fusion_file =  $tumor_filtered_fusion_dir . '/chimeras.filtered.bedpe';
            push @inputs, tumor_unfiltered_fusion_file => $tumor_unfiltered_fusion_file;
            push @inputs, ncbi_human_ensembl_build_id => $ncbi_human_ensembl_build_id;
            push @inputs, tumor_filtered_fusion_file => $tumor_filtered_fusion_file;
            push @dirs, $tumor_filtered_fusion_dir;
            #Check for SV calls file
            if(-e $wgs_build->data_directory . '/effects/svs.hq.annotated'){
                my $wgs_sv_file = $build->wgs_build->data_directory . '/effects/svs.hq.annotated';
                my $tumor_filtered_intersected_fusion_file =  $tumor_filtered_fusion_dir . '/chimeras.filtered.intersected.bedpe';
                push @inputs, wgs_sv_file => $wgs_sv_file;
                push @inputs, tumor_filtered_intersected_fusion_file => $tumor_filtered_intersected_fusion_file;
            }
          }
        }
      }

    #GenerateClonalityPlots
    if ($wgs_build){
      my $clonality_dir = $patient_dir . "/clonality/";
      push @dirs, $clonality_dir;
      push @inputs, clonality_dir => $clonality_dir;      
    }

    #RunCnView
    if ($wgs_build){
      my $cnv_dir = $patient_dir . "/cnv/";
      push @dirs, $cnv_dir;
      push @inputs, cnv_dir => $cnv_dir;
    }

    #SummarizeSvs
    if ($wgs_build){
      my $sv_dir = $patient_dir . "/sv/";
      push @dirs, $sv_dir;
      push @inputs, sv_dir => $sv_dir;
    }

    #CreateMutationSpectrum
    if ($wgs_build) {
        push @inputs, 'wgs_mutation_spectrum_outdir' => $patient_dir . '/mutation-spectrum';
        push @inputs, 'wgs_mutation_spectrum_datatype' => 'wgs';
    }
    if ($exome_build) {
        push @inputs, 'exome_mutation_spectrum_outdir' => $patient_dir . '/mutation-spectrum';
        push @inputs, 'exome_mutation_spectrum_datatype' => 'exome';
    }

    #AnnotateGenesByCategory
    # TODO: most things will use the db object instead of the dir
    # Though it may be better to use the sub-dataset, at least if it can be given a "type" eventually.
    #my $gene_symbol_lists_dir = "/gscmnt/sata132/techd/mgriffit/reference_annotations/GeneSymbolLists/";
    my $gene_symbol_lists_dir = $cancer_annotation_db->data_set_path("GeneSymbolLists");
    push @inputs, 'gene_symbol_lists_dir' => $gene_symbol_lists_dir;
    push @inputs, 'gene_name_columns' => ['mapped_gene_name'];
    push @inputs, 'gene_name_regex' => 'mapped_gene_name';

    # For now it works to create directories here because the data_directory has been allocated.  
    #It is possible that this would not happen until later, which would mess up assigning inputs to many of the commands.
    for my $dir (@dirs) {
        Genome::Sys->create_directory($dir);
    }

    return @inputs;
}

sub _resolve_workflow_for_build {
    # This is called by Genome::Model::Build::start()
    # Returns a Workflow::Operation
    my $self = shift;
    my $build = shift;
    my $lsf_queue = shift;   # TODO: the workflow shouldn't need this yet
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
        igv_session_result
    );

    if ($build->wgs_build) {
        push @output_properties, qw( 
            summarize_wgs_tier1_snv_support_result
            summarize_svs_result
            summarize_cnvs_result
            wgs_mutation_spectrum_result
            clonality_result
            run_cn_view_result
            gene_category_cnv_amp_result
            gene_category_cnv_del_result
            gene_category_cnv_ampdel_result
            gene_category_wgs_snv_result
            gene_category_wgs_indel_result
            dgidb_cnv_amp_result
            dgidb_sv_fusion_result
            dgidb_wgs_snv_result
            dgidb_wgs_indel_result
            wgs_variant_sources_result
        );
    }

    if ($build->exome_build) {
        push @output_properties, qw(
            summarize_exome_tier1_snv_support_result
            exome_mutation_spectrum_result
            gene_category_exome_snv_result
            gene_category_exome_indel_result
            dgidb_exome_snv_result
            dgidb_exome_indel_result
            exome_variant_sources_result
        );
    }

    if ($build->wgs_build and $build->exome_build) {
        push @output_properties, 'summarize_wgs_exome_tier1_snv_support_result';
        push @output_properties, 'gene_category_wgs_exome_indel_result';
        push @output_properties, 'gene_category_wgs_exome_snv_result';
        push @output_properties, 'dgidb_wgs_exome_indel_result';
        push @output_properties, 'dgidb_wgs_exome_snv_result';
    }

    if ($build->wgs_build or $build->exome_build) {
        push @output_properties, 'mutation_diagram_result';
        push @output_properties, 'import_snvs_indels_result';
    }

    if ($build->normal_rnaseq_build){
        push @output_properties, 'normal_tophat_junctions_absolute_result';
        push @output_properties, 'normal_cufflinks_expression_absolute_result';
    }

    if ($build->tumor_rnaseq_build){
        push @output_properties, 'tumor_tophat_junctions_absolute_result';
        push @output_properties, 'tumor_cufflinks_expression_absolute_result';
        push @output_properties, 'gene_category_cufflinks_result';
        push @output_properties, 'gene_category_tophat_result';
        push @output_properties, 'dgidb_cufflinks_result';
        push @output_properties, 'dgidb_tophat_result';
    }
    if ($build->tumor_rnaseq_build){
        if(-e $build->tumor_rnaseq_build->data_directory . '/fusions/chimeras.bedpe'){
            push @output_properties, 'tumor_chimerascan_fusion_filter_result';
            if(-e $build->wgs_build->data_directory . '/effects/svs.hq.annotated'){
                push @output_properties, 'intersect_tumor_fusion_sv_result';
            }
        }
    }
    if ($build->normal_rnaseq_build and $build->tumor_rnaseq_build){
        push @output_properties, 'cufflinks_differential_expression_result';
        push @output_properties, 'gene_category_coding_de_up_result';
        push @output_properties, 'gene_category_coding_de_down_result';
        push @output_properties, 'gene_category_coding_de_result';
    }

    # Make the workflow and some convenience wrappers for adding steps and links

    my $workflow = Workflow::Model->create(
        name => $build->workflow_name,
        input_properties  => \@input_properties, 
        output_properties => \@output_properties,
    );

    my $log_directory = $build->log_directory;
    $workflow->log_dir($log_directory);

    my $input_connector  = $workflow->get_input_connector;
    my $output_connector = $workflow->get_output_connector;

    my %steps_by_name; 
    my $step = 0;
    my $add_step;
    my $add_link;

    $add_step = sub {
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

        # for inputs with these names, always add a link from the input connector
        # should this be default behavior for commands in a model's namespace?
        my $cmeta = $cmd->__meta__;
        for my $pname (qw/
            cancer_annotation_db 
            misc_annotation_db 
            cosmic_annotation_db
        /) {
            my $pmeta = $cmeta->property($pname);
            if ($pmeta) {
                $add_link->($input_connector, $pname, $op, $pname);            
            }
        }
            
        return $op;
    };

    my $converge;

    $add_link = sub {
        my ($from_op, $from_p, $to_op, $to_p) = @_;
        $to_p = $from_p if not defined $to_p;
        
        if (ref($to_p) eq 'ARRAY') {
            Carp::confess("the 'to' property in a link cannot be a list!");
        }

        my $link;
        if (ref($from_p) eq 'ARRAY') {
            my $cname = "Combine: (@$from_p) for \"" . $to_op->name . "\" parameter \'$to_p\'";
            my $converge_op = $converge->($cname,$from_op,$from_p);
            $link = $workflow->add_link(
                left_operation => $converge_op,
                left_property => 'outputs',
                right_operation => $to_op,
                right_property => $to_p
            );
        }
        else {
            $link = $workflow->add_link(
                left_operation => $from_op,
                left_property => $from_p,
                right_operation => $to_op,
                right_property => $to_p
            );
        }
        $link or die "Failed to make link from $link from $from_p to $to_p!";
        return $link;
    };

    my $combo_cnt = 0;
    $converge = sub {
        my $cname = shift;
        my @params1 = @_;
        my @params2 = @_;

        # make a friendly name
        my $name = '';
        my $combo_cnt++;
        my $input_count = 0;
        while (@params1) {
            my $from_op = shift @params1;
            my $from_props = shift @params1;
            unless ($from_props) {
                die "expected \$op2,['p1','p2',...],\$op2,['p3','p4'],...";
            }
            unless (ref($from_props) eq 'ARRAY') {
                die "expected the second param (and every even param) in converge to be an arrayref of property names"; 
            }
            $input_count += scalar(@$from_props);
            if ($name) {
                $name .= " and ";
            }
            else {
                $name = "combine ($combo_cnt): ";
            }
            if ($from_op->name eq 'input_conector') {
                $name = "($combo_cnt)";
            }
            else {
                $name .= $from_op->name . "(";
            }
            for my $p (@$from_props) {
                $name .= $p;
                $name .= "," unless $p eq $from_props->[-1];
            }
            $name .= ")";
        }
        my $op = $workflow->add_operation(
            name => $cname,
            operation_type => Workflow::OperationType::Converge->create(
                input_properties => [ map { "i$_" } (1..$input_count) ],
                output_properties => ['outputs'],
            )
        );
        
        # create links
        my $input_n = 0;
        while (@params2) {
            my $from_op = shift @params2;
            my $from_props = shift @params2;
            for my $from_prop (@$from_props) {
                $input_n++;
                $add_link->(
                    $from_op,
                    $from_prop,
                    $op,
                    "i$input_n"
                );
            }
        }

        # return the op which will have a single "output"
        return $op;
    };
  
    #SummarizeBuilds - Summarize build inputs using SummarizeBuilds.pm
    my $msg = "Creating a summary of input builds using summarize-builds";
    my $summarize_builds_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeBuilds");
    $add_link->($input_connector, 'build_as_array', $summarize_builds_op, 'builds');
    $add_link->($input_connector, 'summarize_builds_outdir', $summarize_builds_op, 'outdir');
    $add_link->($input_connector, 'summarize_builds_skip_lims_reports', $summarize_builds_op, 'skip_lims_reports');
    $add_link->($input_connector, 'summarize_builds_log_file', $summarize_builds_op, 'log_file');
    $add_link->($summarize_builds_op, 'result', $output_connector, 'summarize_builds_result');

    #ImportSnvsIndels - Import SNVs and Indels
    my $import_snvs_indels_op;
    if ($build->wgs_build or $build->exome_build){
      my $msg = "Importing SNVs and Indels from somatic results, parsing, and merging exome/wgs if possible";
      $import_snvs_indels_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::ImportSnvsIndels");
      $add_link->($input_connector, 'import_snvs_indels_outdir', $import_snvs_indels_op, 'outdir');
      $add_link->($input_connector, 'import_snvs_indels_filter_mt', $import_snvs_indels_op, 'filter_mt');
      if ($build->wgs_build){
        $add_link->($input_connector, 'wgs_build', $import_snvs_indels_op, 'wgs_build');
      }
      if ($build->exome_build){
        $add_link->($input_connector, 'exome_build', $import_snvs_indels_op, 'exome_build');
      }
      $add_link->($import_snvs_indels_op, 'result', $output_connector, 'import_snvs_indels_result');
    }

    #GetVariantSources - Determine source variant caller for SNVs and InDels for wgs data
    if ($build->wgs_build) {
      my $msg = "Determining source variant callers of all tier1-3 SNVs and InDels for wgs data";
      my $wgs_variant_sources_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::GetVariantSources");
      $add_link->($input_connector, 'wgs_build_as_array', $wgs_variant_sources_op, 'builds');
      $add_link->($input_connector, 'wgs_variant_sources_dir', $wgs_variant_sources_op, 'outdir');
      $add_link->($wgs_variant_sources_op, 'result', $output_connector, 'wgs_variant_sources_result');
    }
    #GetVariantSources - Determine source variant caller for SNVs and InDels for exome data
    if ($build->exome_build) {
      my $msg = "Determining source variant callers of all tier1-3 SNVs and InDels for exome data";
      my $exome_variant_sources_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::GetVariantSources");
      $add_link->($input_connector, 'exome_build_as_array', $exome_variant_sources_op, 'builds');
      $add_link->($input_connector, 'exome_variant_sources_dir', $exome_variant_sources_op, 'outdir');
      $add_link->($exome_variant_sources_op, 'result', $output_connector, 'exome_variant_sources_result');
    }

    #CreateMutationDiagrams - Create mutation spectrum results for wgs data
    if ($build->wgs_build) {
      my $msg = "Creating mutation spectrum results for wgs snvs using create-mutation-spectrum";
      my $create_mutation_spectrum_wgs_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::CreateMutationSpectrum');
      $add_link->($input_connector, 'wgs_build', $create_mutation_spectrum_wgs_op, 'build');
      $add_link->($input_connector, 'wgs_mutation_spectrum_outdir', $create_mutation_spectrum_wgs_op, 'outdir');
      $add_link->($input_connector, 'wgs_mutation_spectrum_datatype', $create_mutation_spectrum_wgs_op, 'datatype');
      $add_link->($create_mutation_spectrum_wgs_op, 'result', $output_connector, 'wgs_mutation_spectrum_result')
    }

    #CreateMutationDiagrams - Create mutation spectrum results for exome data
    if ($build->exome_build) {
      my $msg = "Creating mutation spectrum results for exome snvs using create-mutation-spectrum";
      my $create_mutation_spectrum_exome_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::CreateMutationSpectrum');
      $add_link->($input_connector, 'exome_build', $create_mutation_spectrum_exome_op, 'build');
      $add_link->($input_connector, 'exome_mutation_spectrum_outdir', $create_mutation_spectrum_exome_op, 'outdir');
      $add_link->($input_connector, 'exome_mutation_spectrum_datatype', $create_mutation_spectrum_exome_op, 'datatype');
      $add_link->($create_mutation_spectrum_exome_op, 'result', $output_connector, 'exome_mutation_spectrum_result')
    }

    #CreateMutationDiagrams - Create mutation diagrams (lolliplots) for all Tier1 SNVs/Indels and compare to COSMIC SNVs/Indels
    if ($build->wgs_build or $build->exome_build) {
      my $msg = "Creating mutation-diagram plots";
      my $mutation_diagram_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::CreateMutationDiagrams");
      if ($build->wgs_build and $build->exome_build) {
          $add_link->($input_connector, ['wgs_build','exome_build'], $mutation_diagram_op, 'builds');
      }
      elsif ($build->wgs_build) {
          $add_link->($input_connector, 'wgs_build', $mutation_diagram_op, 'builds');
      }
      elsif ($build->exome_build) {
          $add_link->($input_connector, 'exome_build', $mutation_diagram_op, 'builds');
      }
      $add_link->($mutation_diagram_op,'result',$output_connector,'mutation_diagram_result');
      
      for my $p (qw/outdir collapse_variants max_snvs_per_file max_indels_per_file/) {
          my $input_name = 'mutation_diagram_' . $p;
          $add_link->($input_connector,$input_name,$mutation_diagram_op,$p);
      }
    }

    #TophatJunctionsAbsolute - Run tophat junctions absolute analysis on normal
    my $normal_tophat_junctions_absolute_op;
    if ($build->normal_rnaseq_build){
      my $msg = "Performing tophat junction expression absolute analysis for normal sample";
      $normal_tophat_junctions_absolute_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::TophatJunctionsAbsolute');
      $add_link->($input_connector, 'normal_rnaseq_build', $normal_tophat_junctions_absolute_op, 'build');
      $add_link->($input_connector, 'normal_tophat_junctions_absolute_dir', $normal_tophat_junctions_absolute_op, 'outdir');
      $add_link->($normal_tophat_junctions_absolute_op, 'result', $output_connector, 'normal_tophat_junctions_absolute_result'); 
    }
    #TophatJunctionsAbsolute - Run tophat junctions absolute analysis on tumor
    my $tumor_tophat_junctions_absolute_op;
    if ($build->tumor_rnaseq_build){
      my $msg = "Performing tophat junction expression absolute analysis for tumor sample";
      $tumor_tophat_junctions_absolute_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::TophatJunctionsAbsolute');
      $add_link->($input_connector, 'tumor_rnaseq_build', $tumor_tophat_junctions_absolute_op, 'build');
      $add_link->($input_connector, 'tumor_tophat_junctions_absolute_dir', $tumor_tophat_junctions_absolute_op, 'outdir');
      $add_link->($tumor_tophat_junctions_absolute_op, 'result', $output_connector, 'tumor_tophat_junctions_absolute_result'); 
    }

    #CufflinksExpressionAbsolute - Run cufflinks expression absolute analysis on normal
    my $normal_cufflinks_expression_absolute_op;
    if ($build->normal_rnaseq_build){
      my $msg = "Performing cufflinks expression absolute analysis for normal sample";
      $normal_cufflinks_expression_absolute_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::CufflinksExpressionAbsolute');
      $add_link->($input_connector, 'normal_rnaseq_build', $normal_cufflinks_expression_absolute_op, 'build');
      $add_link->($input_connector, 'normal_cufflinks_expression_absolute_dir', $normal_cufflinks_expression_absolute_op, 'outdir');
      $add_link->($input_connector, 'cufflinks_percent_cutoff', $normal_cufflinks_expression_absolute_op, 'percent_cutoff');
      $add_link->($normal_cufflinks_expression_absolute_op, 'result', $output_connector, 'normal_cufflinks_expression_absolute_result'); 
    }
    #CufflinksExpressionAbsolute - Run cufflinks expression absolute analysis on tumor
    my $tumor_cufflinks_expression_absolute_op;
    if ($build->tumor_rnaseq_build){
      my $msg = "Performing cufflinks expression absolute analysis for tumor sample";
      $tumor_cufflinks_expression_absolute_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::CufflinksExpressionAbsolute');
      $add_link->($input_connector, 'tumor_rnaseq_build', $tumor_cufflinks_expression_absolute_op, 'build');
      $add_link->($input_connector, 'tumor_cufflinks_expression_absolute_dir', $tumor_cufflinks_expression_absolute_op, 'outdir');
      $add_link->($input_connector, 'cufflinks_percent_cutoff', $tumor_cufflinks_expression_absolute_op, 'percent_cutoff');
      $add_link->($tumor_cufflinks_expression_absolute_op, 'result', $output_connector, 'tumor_cufflinks_expression_absolute_result');
    }

    #CufflinksDifferentialExpression - Run cufflinks differential expression
    my $cufflinks_differential_expression_op;
    if ($build->normal_rnaseq_build and $build->tumor_rnaseq_build){
      my $msg = "Performing cufflinks differential expression analysis for case vs. control (e.g. tumor vs. normal)";
      $cufflinks_differential_expression_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::CufflinksDifferentialExpression');
      $add_link->($input_connector, 'normal_rnaseq_build', $cufflinks_differential_expression_op, 'control_build');
      $add_link->($input_connector, 'tumor_rnaseq_build', $cufflinks_differential_expression_op, 'case_build');
      $add_link->($input_connector, 'cufflinks_differential_expression_dir', $cufflinks_differential_expression_op, 'outdir');
      $add_link->($cufflinks_differential_expression_op, 'result', $output_connector, 'cufflinks_differential_expression_result');
    }

    my $tumor_chimerascan_filtered_fusion_op;
    my $intersect_tumor_fusion_sv_op;
    if ($build->tumor_rnaseq_build){
        #Filter ChimeraScan tumor Fusion output - Use the ChimeraScan::FilterOutput GMT on the tumor fusion calls
        if(-e $build->tumor_rnaseq_build->data_directory . '/fusions/chimeras.bedpe'){
            my $msg = "Filtering tumor ChimeraScan fusion calls.";
            $tumor_chimerascan_filtered_fusion_op = $add_step->($msg, 'Genome::Model::Tools::ChimeraScan::FilterOutput');
            $add_link->($input_connector, 'tumor_unfiltered_fusion_file', $tumor_chimerascan_filtered_fusion_op, 'bedpe_file');
            $add_link->($input_connector, 'ncbi_human_ensembl_build_id', $tumor_chimerascan_filtered_fusion_op, 'annotation_build_id');
            $add_link->($input_connector, 'tumor_filtered_fusion_file', $tumor_chimerascan_filtered_fusion_op, 'output_file');
            $add_link->($tumor_chimerascan_filtered_fusion_op, 'result', $output_connector, 'tumor_chimerascan_fusion_filter_result');
            if(-e $build->wgs_build->data_directory . '/effects/svs.hq.annotated'){
                my $msg = "Intersecting filtered tumor ChimeraScan fusion calls with WGS SV calls.";
                $intersect_tumor_fusion_sv_op = $add_step->($msg, 'Genome::Model::Tools::ChimeraScan::IntersectSv');
                $add_link->($input_connector, 'ncbi_human_ensembl_build_id', $intersect_tumor_fusion_sv_op, 'annotation_build_id');
                $add_link->($input_connector, 'tumor_filtered_intersected_fusion_file', $intersect_tumor_fusion_sv_op, 'output_file');
                $add_link->($input_connector, 'wgs_sv_file', $intersect_tumor_fusion_sv_op, 'sv_output_file');
                $add_link->($tumor_chimerascan_filtered_fusion_op, 'filtered_bedpe_file', $intersect_tumor_fusion_sv_op, 'filtered_bedpe_file');
                $add_link->($intersect_tumor_fusion_sv_op, 'result', $output_connector, 'intersect_tumor_fusion_sv_result');
            }
        }
    }

    #DumpIgvXml - Create IGV xml session files with increasing numbers of tracks and store in a single (WGS and Exome BAM files, RNA-seq BAM files, junctions.bed, SNV bed files, etc.)
    #genome model clin-seq dump-igv-xml --outdir=/gscuser/mgriffit/ --builds=119971814
    $msg = "Create IGV XML session files for varying levels of detail using the input builds";
    my $igv_session_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::DumpIgvXml");
    $add_link->($input_connector, 'build_as_array', $igv_session_op, 'builds');
    $add_link->($input_connector, 'igv_session_dir', $igv_session_op, 'outdir');
    $add_link->($igv_session_op, 'result', $output_connector, 'igv_session_result'); 

    #GenerateClonalityPlots - Run clonality analysis and produce clonality plots
    $msg = "Run clonality analysis and produce clonality plots";
    my $clonality_op;
    if ($build->wgs_build){
      $clonality_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::GenerateClonalityPlots");
      $add_link->($input_connector, 'wgs_build', $clonality_op, 'somatic_var_build');
      $add_link->($input_connector, 'clonality_dir', $clonality_op, 'output_dir');
      $add_link->($input_connector, 'common_name', $clonality_op, 'common_name');
      $add_link->($clonality_op, 'result', $output_connector, 'clonality_result'); 
    }

    #RunCnView - Produce copy number results with run-cn-view.  Relies on clonality step already having been run
    my $run_cn_view_op;
    if ($build->wgs_build){
      $msg = "Use gmt copy-number cn-view to produce copy number tables and images";
      $run_cn_view_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::RunCnView");
      $add_link->($input_connector, 'wgs_build', $run_cn_view_op, 'build');
      $add_link->($input_connector, 'cnv_dir', $run_cn_view_op, 'outdir');
      $add_link->($clonality_op, 'cnv_hmm_file', $run_cn_view_op);
      $add_link->($run_cn_view_op, 'result', $output_connector, 'run_cn_view_result');
    }
   
    #SummarizeCnvs - Generate a summary of CNV results, copy cnvs.hq, cnvs.png, single-bam copy number plot PDF, etc. to the cnv directory
    #This step relies on the generate-clonality-plots step already having been run 
    #It also relies on run-cn-view step having been run already
    if ($build->wgs_build){
        my $msg = "Summarize CNV results from WGS somatic variation";
        my $summarize_cnvs_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeCnvs");
        $add_link->($input_connector, 'cnv_dir', $summarize_cnvs_op, 'outdir');
        $add_link->($input_connector, 'wgs_build', $summarize_cnvs_op, 'build');
        $add_link->($clonality_op, 'cnv_hmm_file', $summarize_cnvs_op);
        $add_link->($run_cn_view_op, 'gene_amp_file', $summarize_cnvs_op);
        $add_link->($run_cn_view_op, 'gene_del_file', $summarize_cnvs_op);
        $add_link->($summarize_cnvs_op, 'result', $output_connector, 'summarize_cnvs_result');
    }

    #SummarizeSvs - Generate a summary of SV results from the WGS SV results
    my $summarize_svs_op;
    if ($build->wgs_build){
        my $msg = "Summarize SV results from WGS somatic variation";
        $summarize_svs_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeSvs");
        $add_link->($input_connector, 'wgs_build_as_array', $summarize_svs_op, 'builds');
        $add_link->($input_connector, 'sv_dir', $summarize_svs_op, 'outdir');
        $add_link->($summarize_svs_op, 'result', $output_connector, 'summarize_svs_result');
    }

    #Add gene category annotations to some output files from steps above. (e.g. determine which SNV affected genes are kinases, ion channels, etc.)
    #Example files to annotate in this way: 
    # AML103/snv/exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv -> gene_category_exome_snv_result
    # AML103/snv/wgs/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv -> gene_category_wgs_snv_result
    # AML103/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv -> gene_category_wgs_exome_snv_result
    # AML103/indel/exome/indels.hq.tier1.v1.annotated.compact.tsv	-> gene_category_exome_indel_result
    # AML103/indel/wgs/indels.hq.tier1.v1.annotated.compact.tsv	-> gene_category_wgs_indel_result
    # AML103/indel/wgs_exome/indels.hq.tier1.v1.annotated.compact.tsv	-> gene_category_wgs_exome_indel_result
    # AML103/cnv/cnview/cnv.All_genes.amp.tsv -> gene_category_cnv_amp_result
    # AML103/cnv/cnview/cnv.All_genes.del.tsv -> gene_category_cnv_del_result
    # AML103/cnv/cnview/cnv.All_genes.ampdel.tsv -> gene_category_cnv_ampdel_result
    # AML103/rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/isoforms.merged.fpkm.expsort.top1percent.tsv -> gene_category_cufflinks_result
    # AML103/rnaseq/tumor/tophat_junctions_absolute/NCBI-human.ensembl-67_37l_v2.Junction.GeneExpression.top1percent.tsv -> gene_category_tophat_result
    # H_MF-F6-8-R12/rnaseq/cufflinks_differential_expression/genes/case_vs_control.coding.hq.up.tsv -> gene_category_coding_de_up_result
    # H_MF-F6-8-R12/rnaseq/cufflinks_differential_expression/genes/case_vs_control.coding.hq.down.tsv -> gene_category_coding_de_down_result
    # H_MF-F6-8-R12/rnaseq/cufflinks_differential_expression/genes/case_vs_control.coding.hq.de.tsv -> gene_category_coding_de_result

    #AnnotateGenesByCategory - gene_category_exome_snv_result
    if ($build->exome_build){
      my $msg = "Add gene category annotations to SNVs identified by exome";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($import_snvs_indels_op, 'exome_snv_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_exome_snv_result');
    }
    #AnnotateGenesByCategory - gene_category_wgs_snv_result
    if ($build->wgs_build){
      my $msg = "Add gene category annotations to SNVs identified by wgs";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($import_snvs_indels_op, 'wgs_snv_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_wgs_snv_result');
    }
    #AnnotateGenesByCategory - gene_category_wgs_exome_snv_result
    if ($build->wgs_build and $build->exome_build){
      my $msg = "Add gene category annotations to SNVs identified by wgs OR exome";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($import_snvs_indels_op, 'wgs_exome_snv_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_wgs_exome_snv_result');
    }
    #AnnotateGenesByCategory - gene_category_exome_indel_result
    if ($build->exome_build){
      my $msg = "Add gene category annotations to InDels identified by exome";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($import_snvs_indels_op, 'exome_indel_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_exome_indel_result');
    }
    #AnnotateGenesByCategory - gene_category_wgs_indel_result
    if ($build->wgs_build){
      my $msg = "Add gene category annotations to InDels identified by wgs";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($import_snvs_indels_op, 'wgs_indel_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_wgs_indel_result');
    }
    #AnnotateGenesByCategory - gene_category_wgs_exome_indel_result
    if ($build->wgs_build and $build->exome_build){
      my $msg = "Add gene category annotations to InDels identified by wgs OR exome";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($import_snvs_indels_op, 'wgs_exome_indel_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_wgs_exome_indel_result');
    }
    #AnnotateGenesByCategory - gene_category_cnv_amp_result
    if ($build->wgs_build){
      my $msg = "Add gene category annotations to CNV amp genes";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($run_cn_view_op, 'gene_amp_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_cnv_amp_result');
    }
    #AnnotateGenesByCategory - gene_category_cnv_del_result
    if ($build->wgs_build){
      my $msg = "Add gene category annotations to CNV del genes";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($run_cn_view_op, 'gene_del_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_cnv_del_result');
    }
    #AnnotateGenesByCategory - gene_category_cnv_ampdel_result
    if ($build->wgs_build){
      my $msg = "Add gene category annotations to CNV amp OR del genes";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($run_cn_view_op, 'gene_ampdel_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_cnv_ampdel_result');
    }
    #AnnotateGenesByCategory - gene_category_cufflinks_result
    if ($build->tumor_rnaseq_build){
      my $msg = "Add gene category annotations to cufflinks top1 percent genes";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($tumor_cufflinks_expression_absolute_op, 'tumor_fpkm_topnpercent_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_cufflinks_result');
    }
    #AnnotateGenesByCategory - gene_category_tophat_result
    if ($build->tumor_rnaseq_build){
      my $msg = "Add gene category annotations to tophat junctions top1 percent genes";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($tumor_tophat_junctions_absolute_op, 'junction_topnpercent_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_tophat_result');
    }
    #AnnotateGenesByCategory - gene_category_coding_de_up_result
    if ($build->normal_rnaseq_build && $build->tumor_rnaseq_build){
      my $msg = "Add gene category annotations to up-regulated coding genes";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($cufflinks_differential_expression_op, 'coding_hq_up_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_coding_de_up_result');
    }
    #AnnotateGenesByCategory - gene_category_coding_de_down_result
    if ($build->normal_rnaseq_build && $build->tumor_rnaseq_build){
      my $msg = "Add gene category annotations to down-regulated coding genes";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($cufflinks_differential_expression_op, 'coding_hq_down_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_coding_de_down_result');
    }
    #AnnotateGenesByCategory - gene_category_coding_de_result
    if ($build->normal_rnaseq_build && $build->tumor_rnaseq_build){
      my $msg = "Add gene category annotations to DE coding genes";
      my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
      $add_link->($cufflinks_differential_expression_op, 'coding_hq_de_file', $annotate_genes_by_category_op, 'infile');
      $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
      $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_coding_de_result');
    }

    #DGIDB gene annotation
    if ($build->exome_build){
        $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'exome_snv_file', $input_connector, $output_connector, 'dgidb_exome_snv_result');
        $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'exome_indel_file', $input_connector, $output_connector, 'dgidb_exome_indel_result');
    }

    if ($build->wgs_build){
        $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'wgs_snv_file', $input_connector, $output_connector, 'dgidb_wgs_snv_result');
        $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'wgs_indel_file', $input_connector, $output_connector, 'dgidb_wgs_indel_result');
        $self->add_dgidb_op_to_flow($add_step, $add_link, $run_cn_view_op, 'gene_amp_file', $input_connector, $output_connector, 'dgidb_cnv_amp_result');
        $self->add_dgidb_op_to_flow($add_step, $add_link, $summarize_svs_op, 'fusion_output_file', $input_connector, $output_connector, 'dgidb_sv_fusion_result');
    }

    if ($build->wgs_build and $build->exome_build) {
        $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'wgs_exome_snv_file', $input_connector, $output_connector, 'dgidb_wgs_exome_snv_result');
        $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'wgs_exome_indel_file', $input_connector, $output_connector, 'dgidb_wgs_exome_indel_result');
    }

    if ($build->tumor_rnaseq_build){
        $self->add_dgidb_op_to_flow($add_step, $add_link, $tumor_cufflinks_expression_absolute_op, 'tumor_fpkm_topnpercent_file', $input_connector, $output_connector, 'dgidb_cufflinks_result');
        $self->add_dgidb_op_to_flow($add_step, $add_link, $tumor_tophat_junctions_absolute_op, 'junction_topnpercent_file', $input_connector, $output_connector, 'dgidb_tophat_result');
    }

    #SummarizeTier1SnvSupport - For each of the following: WGS SNVs, Exome SNVs, and WGS+Exome SNVs, do the following:
    #Get BAM readcounts for WGS (tumor/normal), Exome (tumor/normal), RNAseq (tumor), RNAseq (normal) - as available of course
    #TODO: Break this down to do direct calls to GetBamReadCounts instead of wrapping it.
    for my $run (qw/wgs exome wgs_exome/) {
      if ($run eq 'wgs' and not $build->wgs_build) {
        next;
      }
      if ($run eq 'exome' and not $build->exome_build) {
        next;
      }
      if ($run eq 'wgs_exome' and not ($build->wgs_build and $build->exome_build)) {
        next;
      }
      my $txt_name = $run;
      $txt_name =~ s/_/ plus /g;
      $txt_name =~ s/wgs/WGS/;
      $txt_name =~ s/exome/Exome/;
      $msg = "$txt_name Summarize Tier 1 SNV Support (BAM read counts)";
      my $summarize_tier1_snv_support_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeTier1SnvSupport");
      $add_link->($import_snvs_indels_op, $run . "_snv_file", $summarize_tier1_snv_support_op, $run . "_positions_file");
      $add_link->($input_connector, 'wgs_build', $summarize_tier1_snv_support_op);
      $add_link->($input_connector, 'exome_build', $summarize_tier1_snv_support_op);
      $add_link->($input_connector, 'tumor_rnaseq_build', $summarize_tier1_snv_support_op);
      $add_link->($input_connector, 'normal_rnaseq_build', $summarize_tier1_snv_support_op);
      if ($build->tumor_rnaseq_build){
        $add_link->($tumor_cufflinks_expression_absolute_op, 'tumor_fpkm_file', $summarize_tier1_snv_support_op);
      }
      $add_link->($input_connector, 'verbose', $summarize_tier1_snv_support_op);
      $add_link->($summarize_tier1_snv_support_op, 'result', $output_connector, "summarize_${run}_tier1_snv_support_result");
    }

    # REMINDER:
    # For new steps be sure to add their result to the output connector if they do not feed into another step.
    # When you do that, expand the list of output properties above. 

    my @errors = $workflow->validate();
    if (@errors) {
        for my $error (@errors) {
            $self->error_message($error);
        }
        die "Invalid workflow!";
    }

    return $workflow;
}


sub add_dgidb_op_to_flow {
    my ($self, $add_step, $add_link, $op, $op_prop, $input_connector, $output_connector, $out_prop) = @_;
    my $msg = "Add dgidb gene annotations to $op_prop";
    my $annotate_genes_by_dgidb_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByDgidb");
    $add_link->($op, $op_prop, $annotate_genes_by_dgidb_op, 'input_file');
    $add_link->($input_connector, 'gene_name_regex', $annotate_genes_by_dgidb_op, 'gene_name_regex');
    $add_link->($annotate_genes_by_dgidb_op, 'result', $output_connector, $out_prop);
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
        .*.jpg$
        .*.jpeg$
        .*.png$
        .*/summary/rc_summary.stderr$
        .*._COSMIC.svg$
        .*.clustered.data.tsv$
        .*.SummarizeBuilds.log.tsv$
        .*.DumpIgvXml.log.txt
        .*/mutation_diagrams/cosmic.mutation-diagram.stderr
        .*/mutation_diagrams/somatic.mutation-diagram.stderr
        .*.sequence-context.stderr
        .*.annotate.stderr
        .*/mutation-spectrum/exome/summarize_mutation_spectrum/summarize-mutation-spectrum.stderr
        .*/mutation-spectrum/wgs/summarize_mutation_spectrum/summarize-mutation-spectrum.stderr
        .*LIMS_Sample_Sequence_QC_library.tsv
        .*.R.stderr
    );
};

# Below are some methods to retrieve files from a build
# Ideally they would be tied to creation of these file paths, but they currently 
# throw an exception if the files don't exist. Consider another approach...
sub patient_dir {
    my ($class, $build) = @_;
    unless($build->common_name) {
        die $class->error_message("Common name is not defined.");
    }
    my $patient_dir = $build->data_directory . "/" . $build->common_name;
    unless(-d $patient_dir) {
        die $class->error_message("ClinSeq patient directory not found. Expected: $patient_dir");
    }
    return $patient_dir;
}

sub snv_variant_source_file {
    my ($class, $build, $data_type,) = @_;

    my $patient_dir = $class->patient_dir($build);
    my $source = $patient_dir . "/variant_source_callers";
    my $exception;
    if(-d $source) {
       my $dir = $source . "/$data_type";
       if(-d $dir) {
           my $file = $dir . "/snv_sources.tsv";
           if(-e $file) {
               return $file;
           }
           else {
               $exception = $class->error_message("Expected $file inside $dir and it did not exist.");
           }
       }
       else {
           $exception = $class->error_message("$data_type sub-directory not found in $source."); 
       }
    }
    else {
        $exception = $class->error_message("$source directory not found");
    }
    die $exception;
}

sub clonality_dir {
    my ($class, $build) = @_;

    my $patient_dir = $class->patient_dir($build);
    my $clonality_dir = $patient_dir . "/clonality";
    
    unless(-d $clonality_dir) {
        die $class->error_message("Clonality directory does not exist. Expected: $clonality_dir");
    }
    return $clonality_dir;
}

sub varscan_formatted_readcount_file {
    my ($class, $build) = @_;
    my $clonality_dir = $class->clonality_dir($build);
    my $readcount_file = $clonality_dir ."/allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts.varscan";
    unless(-e $readcount_file) {
        die $class->error_message("Unable to find varscan formatted readcount file. Expected: $readcount_file");
    }
    return $readcount_file;
}

sub cnaseq_hmm_file {
    my ($class, $build) = @_;
    my $clonality_dir = $class->clonality_dir($build);
    my $hmm_file = $clonality_dir . "/cnaseq.cnvhmm";
    unless(-e $hmm_file) {
        die $class->error_message("Unable to find cnaseq hmm file. Expected: $hmm_file");
    }
    return $hmm_file;
}

1;

