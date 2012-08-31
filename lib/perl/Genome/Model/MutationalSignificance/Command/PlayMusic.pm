package Genome::Model::MutationalSignificance::Command::PlayMusic;

use strict;
use warnings;
use Genome;
use Workflow::Simple;

my @has;

my %module_input_exceptions;
my %parallel_by;

my $DONT_USE = "don't use this";

BEGIN {
#This data structure keeps track of parallel_by commands.
    %parallel_by = (
        "Genome::Model::Tools::Music::Bmr::CalcCovgHelper" => "normal_tumor_bam_pair",
        "Genome::Model::MutationalSignificance::Command::CalcCovg" => "somatic_variation_build",
    );
#This data structure specifies the component command modules that are called by
#the workflow, as well as any inputs to those commands that are provided
#inside the workflow rather than as user inputs.  Any input properties from
#the listed modules not in this data structure will be made input properties#to this command.
    %module_input_exceptions = (
        "Genome::Model::Tools::Music::Proximity" => {
        },
        "Genome::Model::Tools::Music::PathScan" => {
            gene_covg_dir => ["Genome::Model::Tools::Music::Bmr::CalcCovg", 'gene_covg_dir'],
            bmr => ["Genome::Model::Tools::Music::Bmr::CalcBmr", 'bmr_output'],
            output_file => ['input connector', 'path_scan_output_file'],
        },
        "Genome::Model::Tools::Music::Smg" => {
            gene_mr_file => ["Genome::Model::Tools::Music::Bmr::CalcBmr", 'gene_mr_file'],
            output_file => ['input connector', 'smg_output_file'],
        },
        "Genome::Model::Tools::Music::MutationRelation" => {
            output_file => ['input connector', 'mutation_relation_output_file'],
            gene_list => ["Genome::Model::Tools::Music::Smg", 'output_file'],
        },
        "Genome::Model::Tools::Music::Bmr::CalcBmr" => {
            output_dir => ["Genome::Model::Tools::Music::Bmr::CalcCovg", 'output_dir'], 
        },
        "Genome::Model::Tools::Music::ClinicalCorrelation" => {
            output_file => ['input connector', 'clinical_correlation_output_file'],
        },
        "Genome::Model::Tools::Music::CosmicOmim" => {
            output_file => ['input connector', 'cosmic_omim_output_file'],
        },
        "Genome::Model::Tools::Music::Pfam" => {
            output_file => ['input connector', 'pfam_output_file'],
        },
        "Genome::Model::Tools::Music::Bmr::CalcCovgHelper" => {
            normal_tumor_bam_pair => ['input connector', 'normal_tumor_bam_pairs'],
            output_dir => ['input connector', 'roi_covg_dir'],
            output_file => $DONT_USE,
        },
        "Genome::Model::Tools::Music::Bmr::CalcCovg" => {
            output_dir => ["Genome::Model::MutationalSignificance::Command::MergeCalcCovg", 'output_dir'],
            cmd_list_file => $DONT_USE,
            cmd_prefix => $DONT_USE, 
            #The following have default values on CalcCovgHelper but not on CalcCovg
            normal_min_depth => ['input connector', 'normal_min_depth'],
            tumor_min_depth => ['input connector', 'tumor_min_depth'],
            min_mapq => ['input connector', 'min_mapq'],
        },
        "Genome::Model::MutationalSignificance::Command::CalcCovg" => {
            somatic_variation_build => ['input connector', 'somatic_variation_builds'],
            output_dir => ['input connector', 'roi_covg_dir'],
            #These are optional on the PlayMusic command, but mandatory on the CalcCovg module
            music_build => ['input connector', 'music_build'], 
            normal_min_depth => ['input connector', 'normal_min_depth'],
            tumor_min_depth => ['input connector', 'tumor_min_depth'],
            min_mapq => ['input connector', 'min_mapq'],
        },
        "Genome::Model::MutationalSignificance::Command::MergeCalcCovg" => {
            output_dir => ['input connector', 'output_dir'],
            output_files => 1, #set source at runtime
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
                unless ($seen{$p->property_name} or exists $module_input_exceptions{$module}{$name}) {
                    my %data = %{ UR::Util::deep_copy($p) };
                    for my $key (keys %data) {
                        delete $data{$key} if $key =~ /^_/;
                    }
                    delete $data{id};
                    delete $data{db_committed};
                    delete $data{class_name};
                    push @has, $name, \%data;
                    $seen{$name} = 1;
                }
            }
        }
    }
}

class Genome::Model::MutationalSignificance::Command::PlayMusic {
    is => ['Command::V2'],
    has => \@has,
    has_optional_input => [
        music_build => {
            is => 'Genome::Model::Build',
            doc => 'Build that is using the software results from this command',
        },
        somatic_variation_builds => {
            is => 'Genome::Model::Build::SomaticVariation',
            is_many => 1,
            doc => 'Builds to analyze.  Each build must match up with one line in bam_list.  Use these for performance enhancement when available',
        },
        log_directory => {
            is => 'Text',
            doc => "Directory to write log files from MuSiC components",
        },
    ],


    has_output => [
        smg_result => {
            is => 'Text',
            doc => 'Output file from Smg tool',
        },
        pathscan_result => {
            is => 'Text',
            doc => 'Output file from PathScan tool',
        },
        mr_result => {
            is => 'Text',
            doc => 'Output file from MutationRelation tool',
        },
        pfam_result => {
            is => 'Text',
            doc => 'Output file from Pfam tool',
        },
        proximity_result => {
            is => 'Text',
            doc => 'Output file from Proximity tool',
        },
        cosmic_result => {
            is => 'Text',
            doc => 'Output file from Cosmic-OMIM tool',
        },
        cct_result => {
            is => 'Text',
            doc => 'Output file from ClinicalCorrelation tool',
        },
    ],
};

sub help_synopsis {
    return <<EOS
This tool takes as parameters all the information required to run the individual tools.  The tools will be run in parallel in an lsf workflow, with values such as bmr passed between tools.  An example usage is:

... genome model mutational-significance play-music\\
        --bam-list input/bams_to_analyze.txt \\
        --numeric-clinical-data-file input/numeric_clinical_data.csv \\
        --maf-file input/myMAF.tsv \\
        --output-dir play_output_dir \\
        --pathway-file input/pathway_db \\
        --reference-sequence input/refseq/all_sequences.fa \\
        --roi-file input/all_coding_regions.bed \\
        --genetic-data-type gene
        --log-directory play_output_dir \\
        --reference-build Build37 \\

EOS
}

sub help_detail {
    return <<EOS
This command can be used to run all of the MuSiC analysis tools on a set of data.  Please see the individual tools for further description of the parameters.
EOS
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    #check somatic_variation_builds and bam_list, if available
    if ($self->somatic_variation_builds) {
        open(BAM_LIST, $self->bam_list);
        my @bam_list;
        while(<BAM_LIST>) {
            push @bam_list, $_;
        }
        close BAM_LIST;
        chomp @bam_list;
        bam_list: foreach my $line (@bam_list) {
            my @fields = split (/\t/, $line);
            my $found = 0;
            foreach my $build ($self->somatic_variation_builds) {
                if ($build->normal_bam eq $fields[1] and $build->tumor_bam eq $fields[2]) {
                    $found = 1;
                    next bam_list;
                }
            }
            if (!$found) {
                $self->error_message("Somatic variation build for bam_list line \"$line\" was not found in the somatic_variation_builds inputs");
                return;
            }
        }
        unless ($self->music_build) {
            $self->error_message("No user build passed in to play-music workflow.  If you pass in somatic variation builds to use as software results, you must also pass in the build that will use the software results");
            return;
        }
    }
    return $self;
}

sub execute {
    my $self = shift;

    my $roi_covg_dir = $self->output_dir."/roi_covgs";
    my $gene_covg_dir = $self->output_dir."/gene_covgs";

    # Create the output directories unless they already exist
    mkdir $roi_covg_dir unless( -e $roi_covg_dir );
    mkdir $gene_covg_dir unless( -e $gene_covg_dir );

    my $workflow = $self->_create_workflow;

    my $meta = $self->__meta__;
    my @all_params = $meta->properties(
        class_name => __PACKAGE__,
        is_input => 1,
    );
    my %properties;
    map {if (defined $self->{$_->property_name}){$properties{$_->property_name} = $self->{$_->property_name}} } @all_params;

    my $bam_list = Genome::Sys->open_file_for_reading($self->bam_list);
    my @normal_tumor_pairs;
    while(my $line = <$bam_list>) {
        chomp $line;
        push @normal_tumor_pairs, $line;
    }

    $properties{normal_tumor_pairs} = \@normal_tumor_pairs;

    $properties{roi_covg_dir} = $roi_covg_dir;

    foreach my $module (qw/path_scan smg mutation_relation clinical_correlation cosmic_omim pfam/) {
        $properties{$module."_output_file"} = $self->output_dir."/$module";
    }

    $workflow->save_to_xml(OutputFile => $self->output_dir . '/play_music_build.xml');

    my @errors = $workflow->validate;
    unless ($workflow->is_valid) {
        die ("Errors encountered while validating workflow\n".join("\n",@errors));
    }

    my $output = Workflow::Simple::run_workflow_lsf(
        $workflow,
        %properties
    );

    unless ( defined $output ) {
        my @errors = @Workflow::Simple::ERROR;
        for (@errors) {
            print STDERR $_->error . "\n";
        }
        return;
    }

    $self->smg_result($output->{smg_result});
    $self->pathscan_result($output->{pathscan_result});
    $self->mr_result($output->{mr_result});
    $self->pfam_result($output->{pfam_result});
    $self->proximity_result($output->{proximity_result});
    $self->cosmic_result($output->{cosmic_result});
    $self->cct_result($output->{cct_result});

    return 1;
}

sub _create_workflow {
    my $self = shift;

    my $meta = $self->__meta__;
    my @input_params = $meta->properties(
        class_name => __PACKAGE__,
        is_input => 1,
    );
    my @input_properties;
    map {if (defined $self->{$_->property_name}){push @input_properties, $_->property_name}} @input_params;

    foreach my $module (qw/path_scan smg mutation_relation clinical_correlation cosmic_omim pfam/) {
        push @input_properties, $module."_output_file";
    }

    push @input_properties, "normal_tumor_pairs";
    push @input_properties, "roi_covg_dir";

    my @output_params = $meta->properties(
        class_name => __PACKAGE__,
        is_output => 1,
    );

    my @output_properties;
    map {push @output_properties, $_->property_name} @output_params;

    my $workflow = Workflow::Model->create(
        name => 'Play Music Inner Workflow',
        input_properties => \@input_properties,
        output_properties =>\@output_properties,
    );

    my $log_directory;
    if ($self->log_directory) {
        $log_directory = $self->log_directory;
    }
    else {
        $log_directory = $self->output_dir;
    }

    $workflow->log_dir($log_directory);

    my $input_connector = $workflow->get_input_connector;
    my $output_connector = $workflow->get_output_connector;

    my $command_module;

    if ($self->somatic_variation_builds) {
        delete $module_input_exceptions{"Genome::Model::Tools::Music::Bmr::CalcCovgHelper"};
        $module_input_exceptions{"Genome::Model::MutationalSignificance::Command::MergeCalcCovg"}{"output_files"} = ['Genome::Model::MutationalSignificance::Command::CalcCovg',"output_file"],
    }
    else {
        delete $module_input_exceptions{"Genome::Model::MutationalSignificance::Command::CalcCovg"};
        $module_input_exceptions{"Genome::Model::MutationalSignificance::Command::MergeCalcCovg"}{"output_files"} = ['Genome::Model::Tools::Music::Bmr::CalcCovgHelper',"final_output_file"],
    }

    my @no_dependencies = ('Genome::Model::Tools::Music::Proximity', 'Genome::Model::Tools::Music::ClinicalCorrelation', 'Genome::Model::Tools::Music::CosmicOmim', 'Genome::Model::Tools::Music::Pfam');
    my @bmr = ('Genome::Model::MutationalSignificance::Command::CalcCovg', 'Genome::Model::Tools::Music::Bmr::CalcCovgHelper', 'Genome::Model::MutationalSignificance::Command::MergeCalcCovg','Genome::Model::Tools::Music::Bmr::CalcCovg', 'Genome::Model::Tools::Music::Bmr::CalcBmr');
    my @depend_on_bmr = ('Genome::Model::Tools::Music::PathScan', 'Genome::Model::Tools::Music::Smg');
    my @depend_on_smg = ('Genome::Model::Tools::Music::MutationRelation');
    for my $command_name (@no_dependencies, @bmr, @depend_on_bmr, @depend_on_smg) {
        $workflow = $self->_append_command_to_workflow($command_name, $workflow)
                    or return;
    }

    my $link;
    $link = $workflow->add_link(
        left_operation => $self->_get_operation_for_module_name('Genome::Model::Tools::Music::Proximity',                                                            $workflow),
        left_property => 'output_file',
        right_operation => $output_connector,
        right_property => 'proximity_result',
    );

    $link = $workflow->add_link(
        left_operation => $self->_get_operation_for_module_name('Genome::Model::Tools::Music::Pfam', $workflow),
        left_property => 'output_file',
        right_operation => $output_connector,
        right_property => 'pfam_result',
    );

    $link = $workflow->add_link(
        left_operation => $self->_get_operation_for_module_name('Genome::Model::Tools::Music::MutationRelation', $workflow),
        left_property => 'output_file',
        right_operation => $output_connector,
        right_property => 'mr_result',
    );

    $link = $workflow->add_link(
        left_operation => $self->_get_operation_for_module_name('Genome::Model::Tools::Music::PathScan', $workflow),
        left_property => 'output_file',
        right_operation => $output_connector,
        right_property => 'pathscan_result',
    );

    $link = $workflow->add_link(
        left_operation => $self->_get_operation_for_module_name('Genome::Model::Tools::Music::Smg', $workflow),
        left_property => 'output_file',
        right_operation => $output_connector,
        right_property => 'smg_result',
    );
    
    $link = $workflow->add_link(
        left_operation => $self->_get_operation_for_module_name('Genome::Model::Tools::Music::CosmicOmim', $workflow),
        left_property => 'output_file',
        right_operation => $output_connector,
        right_property => 'cosmic_result',
    );

    $link = $workflow->add_link(
        left_operation => $self->_get_operation_for_module_name('Genome::Model::Tools::Music::ClinicalCorrelation', $workflow),
        left_property => 'output_file',
        right_operation => $output_connector,
        right_property => 'cct_result',
    );

    my $op = $self->_get_operation_for_module_name('Genome::Model::Tools::Music::Bmr::CalcCovg', $workflow);
    $op->operation_type->lsf_resource("-R \'select[mem>16000] rusage[mem=16000]\' -M 16000000");
    $op = $self->_get_operation_for_module_name('Genome::Model::Tools::Music::Bmr::CalcBmr', $workflow);
    $op->operation_type->lsf_resource("-R \'select[mem>64000] rusage[mem=64000]\' -M 64000000");
    $op = $self->_get_operation_for_module_name('Genome::Model::Tools::Music::Smg', $workflow);
    $op->operation_type->lsf_resource("-R \'select[mem>16000] rusage[mem=16000] span[hosts=1]\' -n ".$self->processors." -M 16000000");
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
    my $command_name = shift;
    my $workflow = shift;
    my $lsf_project = shift;
    my $lsf_queue = shift;
    my $command_module = $command_name;
    unless (defined $module_input_exceptions{$command_name}) {
        return $workflow;
    }
    my $command_meta = $command_module->__meta__;
    my $operation_name = $command_module;
    my $operation;
    if ($parallel_by{$operation_name}){
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
        if(!$property->is_optional or defined $property_def) {
            if (!$property->is_optional and not defined $property_def) {
                die ("Non-optional property ".$property->property_name." is not provided\n");
            }
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

        }
    }
    return $workflow;
}

1;

