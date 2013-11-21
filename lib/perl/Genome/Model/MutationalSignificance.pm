package Genome::Model::MutationalSignificance;

use strict;
use warnings;

use Genome;
use Genome::Model::ClinSeq::Util;

my @has_param;
my %module_input_exceptions;
my %parallel_by;

my $DONT_USE;

BEGIN {
    $DONT_USE = "don't use this";
    %parallel_by = (
        "Genome::Model::MutationalSignificance::Command::CreateMafFile" => "somatic_variation_build",
    );
    %module_input_exceptions = (
        'Genome::Model::MutationalSignificance::Command::MergeMafFiles' => {
            maf_files => ["Genome::Model::MutationalSignificance::Command::CreateMafFile", 'maf_file'],
            maf_path => ["input connector", "merged_maf_path"],
        },
        'Genome::Model::MutationalSignificance::Command::CreateMafFile' => {
            somatic_variation_build => ['input connector', 'somatic_variation_builds'],
            output_dir => ['input connector', 'create_maf_output_dir'],
            cosmic_dir => ["input connector", 'cosmic_dir'],
            review_file_dir => ["input connector", 'review_file_dir'],
            regulatory_columns_to_check => ["input connector", "regulatory_columns_to_check"],
        },
        'Genome::Model::MutationalSignificance::Command::CreateROI' => {
            annotation_build => ['input connector', 'annotation_build'], #input to model, not param
            extra_rois => ['input connector', 'extra_rois'], #TODO: get from somatic variation models?
            regulome_bed => ['input connector', 'regulome_bed'],
            include_ensembl_annot => ['input connector', 'include_ensembl_annot'],
        },
        'Genome::Model::MutationalSignificance::Command::CreateBamList' => {
            somatic_variation_builds => ['input connector', 'somatic_variation_builds'],
            bam_list => ["input connector", "bam_list"],
        },
        'Genome::Model::MutationalSignificance::Command::PlayMusic' => {
            input_clinical_correlation_matrix_file => $DONT_USE,
            bam_list => ["Genome::Model::MutationalSignificance::Command::CreateBamList", "bam_list"],
            maf_file => ["Genome::Model::MutationalSignificance::Command::MergeMafFiles", "maf_path"],
            roi_file => ["Genome::Model::MutationalSignificance::Command::CreateROI", "roi_path"],
            reference_sequence => ['input connector', 'reference_sequence'],
            output_dir => ['input connector', 'output_dir'],
            log_directory => ['input connector', 'log_directory'],
            music_build => ['input connector', 'music_build'],
            somatic_variation_builds => ['input connector', 'somatic_variation_builds'],
            pathway_file => ['input connector', 'pathway_file'],
            categorical_clinical_data_file => ['input connector', 'categorical_clinical_data_file'],
            numeric_clinical_data_file => ['input connector', 'numeric_clinical_data_file'],
            glm_clinical_data_file => ['input connector', 'glm_clinical_data_file'],
            glm_model_file => ['input connector', 'glm_model_file'],
            omimaa_dir => ['input connector', 'omimaa_dir'],
            cosmic_dir => ['input connector', 'cosmic_dir'],
            bmr_modifier_file => ['input connector', 'bmr_modifier_file'],
            mutation_matrix_file => ['input connector', 'mutation_matrix_file'],
            clinical_correlation_matrix_file => ['input connector', 'clinical_correlation_matrix_file'],
            reference_build => ['input connector', 'reference_build'],
        },
        'Genome::Model::MutationalSignificance::Command::RunReports' => {
            maf_file => ["Genome::Model::MutationalSignificance::Command::MergeMafFiles", "maf_path"],
            output_dir => ['input connector', 'output_dir'],
            annotation_build => ['input connector', 'annotation_build'],
            somatic_variation_builds => ['input connector', 'somatic_variation_builds'],
        },
        'Genome::Model::MutationalSignificance::Command::CompileValidationList' => {
            significantly_mutated_gene_list => $DONT_USE,
            use_tier_1 => ['input connector', "use_tier_1"],
            use_tier_2 => ['input connector', "use_tier_2"],
            use_tier_3 => ['input connector', "use_tier_3"],
            use_tier_4 => ['input connector', "use_tier_4"],
            significant_variant_list => ['input connector', 'significant_variant_list'],
            somatic_variation_builds => ['input connector', 'somatic_variation_builds'],
            reference_sequence_build => ['input connector', 'reference_sequence_build'],
            exon_bed => ["Genome::Model::MutationalSignificance::Command::CreateROI", "roi_path"],
            regions_of_interest => ['input connector', 'regions_of_interest'],
            gene_black_lists => ['input connector', 'gene_black_lists'],
            additional_snv_lists => ['input connector', 'additional_snv_lists'],
            additional_indel_lists => ['input connector', 'additional_indel_lists'],
            additional_sv_lists => ['input connector', 'additional_sv_lists'],
            variant_black_lists => ['input connector', 'variant_black_lists'],
        },
    );
    my %additional_params = (
        play_music => {
            is => 'Boolean',
            default => 0,
            doc => "Whether to run Play Music as a part of determining which mutations are significant",
        },
    );
    my %seen;
    my @modules = keys %module_input_exceptions;
    foreach my $module (@modules) {
        my $module_meta = UR::Object::Type->get($module);
        my @p = $module_meta->properties;
        for my $p (@p) {
            if ($p->can("is_input") and $p->is_input){
                my $name = $p->property_name;
                unless ($seen{$name} or $module_input_exceptions{$module}{$name}) {
                    my %data = %{$p};
                    for my $key (keys %data) {
                        delete $data{$key} if $key =~ /^_/;
                    }
                    delete $data{id};
                    delete $data{db_committed};
                    delete $data{class_name};
                    delete $data{is_input};
                    $data{is_param} = 1;
                    push @has_param, $name, \%data;
                    $seen{$name} = 1;
                }
            }
        }    
    }
    foreach my $param (keys %additional_params) {
        push @has_param, $param, $additional_params{$param};
    }
}

class Genome::Model::MutationalSignificance {
    is        => 'Genome::Model',
    has => [
        processing_profile => { is => 'Genome::ProcessingProfile', id_by => 'processing_profile_id', },
    ],
    has_input => [
        clinseq_models => {
            is    => 'Genome::Model::ClinSeq',
            is_many => 1,
            doc => 'clinseq models to evaluate',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'annotation to use for roi file',
        },
        review_file_dir => {
            is => 'String',
            is_optional => 1,
            doc => 'Path to directory of variant files with reviews.  Any variant with a review status other than S or V will be  ignored.',
        },
        regions_of_interest => {
            is => 'Genome::FeatureList',
            is_many => 1,
            doc => 'FeatureLists of regions to include in validation',
            is_optional => 1,
        },
        gene_black_lists => {
            is => 'File',
            is_many => 1,
            is_optional => 1,
            doc => 'Lists of genes to exclude from the validation.  One gene symbol per line.  Gene symbols must match symbols in annotation files',
        },
        variant_black_lists => {
            is => 'Genome::FeatureList',
            is_many => 1,
            is_optional => 1,
            doc => 'FeatureLists of variants to exclude from the validation.',
        },
        additional_snv_lists => {
            is => 'Genome::FeatureList',
            is_many => 1,
            is_optional => 1,
            doc => 'FeatureLists of additional snv variants to include in validation',
        },
        additional_indel_lists => {
            is => 'Genome::FeatureList',
            is_many => 1,
            is_optional => 1,
            doc => 'FeatureLists of additional indel variants to include in validation',
        },
        additional_sv_lists => {
            is => 'Genome::FeatureList',
            is_many => 1,
            is_optional => 1,
            doc => 'FeatureLists of additional structural variants to include in validation',
        },
        pathway_file => {
            is => 'File',
            is_optional => 1,
            doc => "Tab-delimited file of pathway information",
            default_value => "/gscmnt/gc2108/info/medseq/ckandoth/music/pathscan/all_pathway_files/KEGG_120910",
        },
        categorical_clinical_data_file => {
            is => 'File',
            is_optional => 1,
            doc => "Table of samples (y) vs. categorical clinical data category (x)",
        },
        numeric_clinical_data_file => {
            is => 'File',
            is_optional => 1,
            doc => 'Table of samples (y) vs. numeric clinical data category (x)',
        },
        glm_clinical_data_file => {
            is => 'File',
            is_optional => 1,
            doc => "Clinical traits, mutational profiles, other mixed clinical data",
        },
        glm_model_file => {
            is => 'File',
            is_optional => 1,
            doc => "File outlining the type of model, response variable, covariants, etc. for the GLM analysis.",
        },
        omimaa_dir => {
            is => 'File',
            is_optional => 1,
            doc => "omim amino acid mutation database folder",
            default_value => "/gsc/scripts/opt/genome/db/omim/20110721",
        },
        cosmic_dir => {
            is => 'File',
            is_optional => 1,
            doc => "cosmic amino acid mutation database folder",
            default_value => "/gsc/scripts/opt/genome/db/cosmic/54",
        },
        bmr_modifier_file => {
            is => 'File',
            is_optional => 1,
            doc => "Tab delimited multipliers per gene that modify BMR before testing",
        },
        extra_rois => {
            is => 'Genome::FeatureList',
            is_optional => 1,
            is_many => 1,
            doc => 'Extra ROI files to include in MuSiC analysis',
        },
        regulome_bed => {
            is => 'Genome::FeatureList',
            is_optional => 1,
            doc => 'Bed file with regulome db scores for regions',
        },
    ],
    has_param => \@has_param,
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
TODO
EOS
}

sub help_detail_for_create_profile {
    return <<EOS
genome processing-profile create mutational-significance --based-on "Oct 2012 Default" --name "New example processing profile with larger flank size" --flank-size 2
EOS
}

sub _help_detail_for_model_define {
    return <<EOS
The mutational significance genome model represents analysis on a collection of clinseq models.  There is an option to run the MuSiC pipeline to run various statistical tests on the collection of variants.  Additionally, a list of variants can be compiled that should be suitable for use as a validation list.
EOS
}

sub _help_synopsis_for_model_define {
    return <<EOS
First, determine what group of ClinSeq models you want to work on.  This might be all models in a model group, for example.
Next, create a population group to use as subject:
genome population-group create --models model_groups.id=29624
Finally, define your model using the query for ClinSeq models and the id of the population group you just created.
genome model define mutational-significance --model-name "Example Mutational Significance" --processing-profile "Oct 2012 Default" --annotation-build 106409619 --clinseq-models model_groups.id=29624 --subject id=2882703149
EOS
}

sub _resolve_workflow_for_build {

    # This is called by Genome::Model::Build::start()
    # Returns a Workflow::Operation
    # By default, builds this from stages(), but can be overridden for custom workflow.
    my $self = shift;
    my $build = shift;
    my $lsf_queue = shift; # TODO: the workflow shouldn't need this yet
    my $lsf_project = shift;
     
    my @input_properties;    

    if (!defined $lsf_queue || $lsf_queue eq '' || $lsf_queue eq 'inline') {
        $lsf_queue = $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT};
    }
    if (!defined $lsf_project || $lsf_project eq '') {
        $lsf_project = 'build' . $build->id;
    }
    my $meta = $self->__meta__;
    my $build_meta = $build->__meta__;
    my @input_params = $build_meta->properties(        
        class_name => $build->class,        
        is_input => 1,
    );

    map {my $name = $_->property_name; my @values = $build->$name; if (defined $values[0]){push @input_properties,   $name} }@input_params;

    @input_params = $meta->properties(
        class_name => $self->class,
        is_param => 1,
    );

    map {my $name = $_->property_name; my @values =  $self->$name; if (defined $values[0]){push @input_properties,   $name} }@input_params;

    push @input_properties, "output_dir";
    push @input_properties, "music_build";
    push @input_properties, "clinical_data_file";
    push @input_properties, "merged_maf_path";
    push @input_properties, "create_maf_output_dir";
    push @input_properties, "bam_list";
    push @input_properties, "reference_sequence";
    push @input_properties, "reference_sequence_build";
    push @input_properties, "log_directory";
    push @input_properties, "significant_variant_list";
    push @input_properties, "mutation_matrix_file";
    #push @input_properties, "categorical_clinical_data_file";
    #push @input_properties, "numeric_clinical_data_file";
    push @input_properties, "clinical_correlation_matrix_file";
    push @input_properties, "reference_build";
    push @input_properties, "regulatory_columns_to_check";
    push @input_properties, "somatic_variation_builds";

    my @output_properties;
    if ($self->play_music) {
        if ($self->run_smg or $self->run_mutation_relation) {
            push @output_properties, 'smg_result';
        }
        if ($self->run_path_scan) {
            push @output_properties, 'pathscan_result';
        }
        if ($self->run_mutation_relation) {
            push @output_properties, 'mr_result';
        }
        if ($self->run_pfam) {
            push @output_properties, 'pfam_result';
        }
        if ($self->run_proximity) {
            push @output_properties, 'proximity_result';
        }
        if ($self->run_cosmic_omim) {
            push @output_properties, 'cosmic_result';
        }
        if ($self->run_clinical_correlation) {
            push @output_properties, 'cct_result';
        }
        if ($self->run_reports) {
            push @output_properties, 'output_dir';
        }
    }
    else {
        @output_properties = ('roi_path', 'maf_path', 'bam_list');
    }
    push @output_properties, "significant_variant_list";

    my $workflow = Workflow::Model->create(
        name => $build->workflow_name,
        input_properties => \@input_properties,
        output_properties => \@output_properties,
    );

    my $log_directory = $build->log_directory;
    $workflow->log_dir($log_directory);
 
    my $output_connector = $workflow->get_output_connector;

    my @commands = ('Genome::Model::MutationalSignificance::Command::CreateMafFile','Genome::Model::MutationalSignificance::Command::MergeMafFiles','Genome::Model::MutationalSignificance::Command::CreateROI','Genome::Model::MutationalSignificance::Command::CreateBamList','Genome::Model::MutationalSignificance::Command::CompileValidationList');
    if ($self->run_reports) {
        push @commands, 'Genome::Model::MutationalSignificance::Command::RunReports';
    }

    for my $command_name (@commands) {
        $workflow = $self->_append_command_to_workflow($command_name,
            $workflow, $lsf_project, $lsf_queue) or return;
    }
    my $link;
    if ($self->play_music) {

        $workflow = $self->_append_command_to_workflow("Genome::Model::MutationalSignificance::Command::PlayMusic", $workflow, $lsf_project, $lsf_queue) or return;
        if ($self->run_proximity) {
            $link = $workflow->add_link(
                left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::PlayMusic',
                    $workflow),
                left_property => 'proximity_result',
                right_operation => $output_connector,
                right_property => 'proximity_result',
            );
        }
        if ($self->run_pfam) {
            $link = $workflow->add_link(
                left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::PlayMusic',
                    $workflow),
                left_property => 'pfam_result',
                right_operation => $output_connector,
                right_property => 'pfam_result',
            );
        }
        if ($self->run_mutation_relation) {
            $link = $workflow->add_link(
                left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::PlayMusic',
                    $workflow),
                left_property => 'mr_result',
                right_operation => $output_connector,
                right_property => 'mr_result',
            );
        }
        if ($self->run_path_scan) {
            $link = $workflow->add_link(
                left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::PlayMusic',
                    $workflow),
                left_property => 'pathscan_result',
                right_operation => $output_connector,
                right_property => 'pathscan_result',
            );
        }
        if ($self->run_smg or $self->run_mutation_relation) {
            $link = $workflow->add_link(
                left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::PlayMusic',
                    $workflow),
                left_property => 'smg_result',
                right_operation => $output_connector,
                right_property => 'smg_result',
            );
        }
        if ($self->run_cosmic_omim) {
            $link = $workflow->add_link(
                left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::PlayMusic',
                    $workflow),
                left_property => 'cosmic_result',
                right_operation => $output_connector,
                right_property => 'cosmic_result',
            );
        }

        if ($self->run_clinical_correlation) {
            $link = $workflow->add_link(
                left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::PlayMusic',
                    $workflow),
                left_property => 'cct_result',
                right_operation => $output_connector,
                right_property => 'cct_result',
            );
        }
        if ($self->run_smg or $self->run_mutation_relation) {
            $link = $workflow->add_link(
                left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::PlayMusic',
                    $workflow),
                left_property => 'smg_result',
                right_operation => $self->_get_operation_for_module_name("Genome::Model::MutationalSignificance::Command::CompileValidationList",
                    $workflow),
                right_property => "significantly_mutated_gene_list",
            );
        }
        if ($self->run_reports) {
            $link = $workflow->add_link(
                left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::PlayMusic', $workflow),
                left_property => 'output_dir',
                right_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::RunReports', $workflow),
                right_property => 'output_dir',
            );
            $link = $workflow->add_link(
                left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::RunReports', $workflow),
                left_property => 'output_dir',
                right_operation => $output_connector,
                right_property => 'output_dir',
            );
        }
    }
    else {
        $link = $workflow->add_link(
            left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::CreateBamList', $workflow),
            left_property => 'bam_list',
            right_operation => $output_connector,
            right_property => 'bam_list',
        );
        $link = $workflow->add_link(
            left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::MergeMafFiles', $workflow),
            left_property => 'maf_path',
            right_operation => $output_connector,
            right_property => 'maf_path',
        );
        $link = $workflow->add_link(
            left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::CreateROI', $workflow),
            left_property => 'roi_path',
            right_operation => $output_connector,
            right_property => 'roi_path',
        );
    }
    $link = $workflow->add_link(
        left_operation => $self->_get_operation_for_module_name('Genome::Model::MutationalSignificance::Command::CompileValidationList', $workflow),
        left_property => 'significant_variant_list',
        right_operation => $output_connector,
        right_property => 'significant_variant_list',
    );

    my @errors = $workflow->validate;
    
    unless ($workflow->is_valid) {
        my $message = join("\n", "Workflow is not valid:", @errors);
        $self->error_message($message);
        return;
    }

    return $workflow;
}

sub _get_operation_for_module_name {
    my $self = shift;
    my $operation_name = shift;
    my $workflow = shift;

    foreach my $op ($workflow->operations) {
        if ($op->name eq $operation_name) {
            return $op;
        }
    }
    return;
}

sub _append_command_to_workflow {
    my $self = shift;
    my $command_module = shift;
    my $workflow = shift;
    my $lsf_project = shift;
    my $lsf_queue = shift;
    my $command_meta = $command_module->__meta__;
    my $operation_name = $command_module;
    my $operation;
    if ($parallel_by{$operation_name}) {
        $operation = $workflow->add_operation(
            name => $operation_name,
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => $command_module,
            ),
            parallel_by => $parallel_by{$operation_name},
        );
    }
    else {
        $operation = $workflow->add_operation(
            name => $operation_name,
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => $command_module,
            )
        );
    }
    $operation->operation_type->lsf_queue($lsf_queue);
    $operation->operation_type->lsf_project($lsf_project);
    for my $property ($command_meta->_legacy_properties()) {
        next unless exists $property->{is_input} and $property->{is_input};
        my $property_name = $property->property_name;
        my $property_def = $module_input_exceptions{$operation_name}{$property_name};
        if (defined $property_def and $property_def eq $DONT_USE) {
            $property_def = undef;
        }
        if (!$property_def) {
            if (grep {/^$property_name$/} @{$workflow->operation_type->input_properties}) {
                $property_def = [$workflow->get_input_connector->name, $property_name];
            }
        }
        if(!$property->is_optional) {
            if (not defined $property_def) {
                die ("Non-optional property ".$property->property_name." is not provided\n");
            }
            $workflow = $self->_add_link($property_name, $property_def, $operation, $workflow);
        }
        elsif (defined $property_def) { 
            if ($property_def->[0] eq $workflow->get_input_connector->name) {
                if (grep {/^$property_name$/} @{$workflow->operation_type->input_properties}) {
                    $workflow = $self->_add_link($property_name, $property_def, $operation, $workflow);
                }
            }
            else {
                $workflow = $self->_add_link($property_name, $property_def, $operation, $workflow);
            }
        }
    }
    return $workflow;
}

sub _add_link {
    my $self = shift;
    my $property_name = shift;
    my $property_def = shift;
    my $operation = shift;
    my $workflow = shift;

    my $from_op = $self->_get_operation_for_module_name($property_def->[0], $workflow);
    if (!$from_op) {
        print "looking for left operation ".$property_def->[0]."\n";
        print "left property ".$property_def->[1]."\n";
        print "right operation ".$operation->name."\n";
        print "right property ".$property_name."\n";
        die ("Didn't get a from operation for the link\n");
    }
    my $link = $workflow->add_link(
        left_operation => $from_op,
        left_property => $property_def->[1],
        right_operation => $operation,
        right_property => $property_name,
    );
    return $workflow;
}

sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    my %inputs;
 
    my $meta = $self->__meta__;
    my $build_meta = $build->__meta__;
    my @all_params = $build_meta->properties(
        class_name => $build->class,
        is_input => 1,
    );
    
    map {my $name = $_->property_name; my @values = $build->$name; if (defined $values[0]){my $value; if ($_->is_many){$value = \@values} else{$value = $values[0]} $inputs{$name} = $value}} @all_params;
    
    @all_params = $meta->properties(
        class_name => $self->class,
        is_param => 1,
    );

    map {my $name = $_->property_name; my @values = $self->$name; if (defined $values[0]){my $value; if ($_->is_many){$value = \@values} else{$value = $values[0]} $inputs{$name} = $value}} @all_params;
 
    my @builds = $build->clinseq_builds;
    my $base_dir = $build->data_directory;

    if ($build->review_file_dir) {
        $inputs{review_file_dir} = $build->review_file_dir;
    }
    my @extra_rois = $build->extra_rois;
    if (@extra_rois) {
        $inputs{extra_rois} = [@extra_rois];
        $inputs{regulatory_columns_to_check} = [map {$_->name} @extra_rois];
    }
    if ($build->regulome_bed) {
        $inputs{regulome_bed} = $build->regulome_bed;
    }
    $inputs{music_build} = $build;
    $inputs{log_directory} = $build->log_directory;
    $inputs{merged_maf_path} = $base_dir."/final.maf";
    $inputs{create_maf_output_dir} = $base_dir;
    $inputs{reference_sequence_build} = Genome::Model::ClinSeq::Util::resolve_reference_sequence_build($builds[0]);
    $inputs{reference_sequence} = $inputs{reference_sequence_build}->fasta_file;
    $inputs{output_dir} = $base_dir;
    $inputs{bam_list} = $base_dir."/bam_list.txt";
    $inputs{clinical_data_file} = $base_dir."/clinical_data.txt";
    $inputs{significant_variant_list} = $base_dir."/significant_variant_list";
    $inputs{mutation_matrix_file} = $base_dir."/mutation_matrix_file";
    $inputs{clinical_correlation_matrix_file} = $base_dir."/clinical_correlation_matrix_file";
    $inputs{somatic_variation_builds} = [map {$_->wgs_build} @builds]; #TODO:get exome build as well?
    #$inputs{categorical_clinical_data_file} = $build->categorical_clinical_data_file;
    #$inputs{numeric_clinical_data_file} = $build->numeric_clinical_data_file;

    my $reference_build_name = $inputs{reference_sequence_build}->name;
    my $calculated_reference_build = "37";
    unless ($reference_build_name =~ /$calculated_reference_build/) {
        $calculated_reference_build = "36";
    }

    $inputs{reference_build} = "Build".$calculated_reference_build;

    return %inputs;
}

1;
