package Genome::Model::Tools::Relationship::BackfillPolymuttVcf;

use strict;
use warnings;
use Genome;
use File::Basename;
use Workflow::Simple;

class Genome::Model::Tools::Relationship::BackfillPolymuttVcf {
    is => 'Genome::Command::Base',
    doc => 'Use polymutt to force-genotype all segregating sites and produce one final vcf',
    has_input => [
       model_group => {
            is => 'Genome::ModelGroup',
       },
       output_dir=> {
           is => "Text",
           doc=>"place to store subsetted glfs and their resulting vcfs."
       },
       chr2process=> {
           doc => "Chromsomes to process. Currently not optional since polymutt will rarely work without it (due to GL contig mismatching",
       },
       use_qc_failed_builds => {
            doc => "If set to true, ignore QC failed status in builds when gathering builds from models",
            default => 1,
       }
    ],
    has_optional_input => [
       segregating_sites_file => {
           is => "Text",
           doc => "Optional roi file (just chrom pos) to use instead of generating one from segregating sites in the model group",
       },
       polymutt_version=> {
           is => "Text",
           default => '0.13',
           doc => "Polymutt version to use when force genotyping. Needs to have the --pos option available, which v0.11 does not",
       },
       joinx_version=> {
           is => "Text",
           default => '1.7',
       },
       # Roi limiting params
       roi_file => {
            is => 'Text',
            doc => 'Set this along with roi_name to limit the original polymutt vcfs (from the model group) to roi target regions',
        },
        roi_name => {
            is => 'Text',
            doc => 'Set this along with roi_file to limit the original polymutt vcfs (from the model group) to roi target regions',
        },
        wingspan => {
            is => 'Text',
            doc => 'Set this to add a wingspan to region limiting',
        },
    ],
    has_transient_optional => [
        _builds => {
            is_many => 1,
            is => "Genome::Model::Build::PhenotypeCorrelation",
        },
    ],
    has_param => [
        lsf_resource => {
            is => 'Text',
            default => "-R 'span[hosts=1] rusage[mem=1000] -n 4'",
        },
        lsf_queue => {
            is => 'Text',
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
    ],
};

sub help_detail {
    'Use polymutt to force-genotype all segregating sites and produce one final vcf'
}

#/gscuser/dlarson/src/polymutt.0.01/bin/polymutt -p 20000492.ped -d 20000492.dat -g 20000492.glfindex --minMapQuality 1 --nthreads 4 --vcf 20000492.standard.vcf
sub execute {
    my $self=shift;

    $self->validate_inputs;

    $self->status_message("Resolving valid builds for the given model group");
    $self->resolve_valid_builds;

    $self->status_message("Preparing directories");
    $self->prepare_directories;

    $self->status_message("Performing region limiting on the original polymutt variant vcfs if requested");
    $self->region_limit_model_group;

    $self->status_message("Resolving the segregating sites file to be used");
    my $segregating_sites_file = $self->resolve_segregating_sites_file;

    $self->status_message("Running Polymutt with the segregating sites file");
    $self->run_polymutt_on_segregating_sites($segregating_sites_file);


    $self->status_message("Backfilling all original vcfs with force-genotype data");
    $self->combine_individual_vcfs;

    $self->status_message("Combining all backfilled vcfs into one final vcf");
    $self->create_final_vcf;

    return 1;
}

sub validate_inputs {
    my $self = shift;

    if ($self->roi_file || $self->roi_name || $self->wingspan) {
        unless ($self->roi_file && $self->roi_name&& $self->wingspan) {
            die $self->error_message("roi_file, roi_name and wingspan must all be defined together if one of them is provided");
        }
    }

    return 1;
}

sub resolve_valid_builds {
    my $self = shift;
    my $mg = $self->model_group;
    my @builds;

    my $no_builds = 0;
    my $bad_builds = 0;
    my $ok = 0;
    my @models = $mg->models;
    for my $model (@models) {
        my $build = $model->last_succeeded_build;
        if($build) {
            if ($build->qc_succeeded || $self->use_qc_failed_builds) {
                push @builds, $build;
                $ok++;
            } else {
                $self->status_message("Skipping " . $model->name . " (" . $model->id  . ") because builds failed QC and use_qc_failed_builds is not set");
                $bad_builds++;
            }
        } else {
            $self->status_message("Skipping " . $model->name . " (" . $model->id  . ") because there are no successful builds.");
            $no_builds++;
        }
    }

    $self->status_message("Out of " . scalar(@models) . ", $ok are fine, $bad_builds were skipped because of QC, and $no_builds were skipped due to lack of succeeded builds");

    $self->_builds(\@builds);

    return 1;
}

# Create the output dir and all family subdirs. Prepare the family subdirs by symlinking the glf files and glf index.
sub prepare_directories {
    my $self = shift;

    Genome::Sys->create_directory($self->output_dir);
    my @builds = $self->_builds;
    for my $build (@builds) {
        my $subdir = $self->subdir_for_build($build);
        Genome::Sys->create_directory($subdir);

        my @input_glfs = $self->glfs_for_build($build);
        my $ped_file = $build->pedigree_file_path->path;
        my $glf_index = $self->new_glf_index_for_build($build);

        my @output_glfs;
        for my $input_glf (@input_glfs) {
            my $output_glf = "$subdir/" . basename($input_glf);
            Genome::Sys->create_symlink($input_glf, $output_glf);
            push @output_glfs, $output_glf;
        }

        Genome::Sys->create_symlink($self->old_glf_index_for_build($build), $self->new_glf_index_for_build($build));
    }

    return 1;
}

sub resolve_segregating_sites_file {
    my $self = shift;

    my $segregating_sites_file = $self->segregating_sites_file;
    unless($segregating_sites_file) {
        $self->status_message("No segregating sites file supplied, attempting to make a segregating sites file of all possible sites from Model group outputs...");
        $segregating_sites_file = $self->assemble_list_of_segregating_sites();
    }
    unless (-s $segregating_sites_file) {
        die $self->error_message("segregating sites file $segregating_sites_file does not exist or is empty");
    }

    return $segregating_sites_file;
}

# Gather a list of all variant sites from all builds
# FIXME can probably be replaced with a call to joinx vcf-merge piped to cut -f1,2 to produce chr pos
sub assemble_list_of_segregating_sites {
    my $self = shift;

    my $output_dir = $self->output_dir;
    my @builds = $self->_builds;
    my $segregating_sites_filename  = $output_dir . "/segregating_sites.bed";

    my %positions;
    #open every final original snvs.vcf.gz and accumulate all unique positions that will need to be polled across the model group. (CHR POS)
    for my $build (@builds) {
        my $snvs_vcf = $self->variant_vcf_for_build($build);
        my $fh = Genome::Sys->open_gzip_file_for_reading($snvs_vcf);
        while(my $line = $fh->getline) {
            if($line =~ m/^#/) {
                next;
            } else {
                my ($chr, $pos, @rest) = split "\t", $line;
                $positions{$chr}{$pos} = 1;
            }
        }
        $fh->close();
    }

    #print segregating sites file
    my $sites_output_fh = Genome::Sys->open_file_for_writing($segregating_sites_filename);
    for my $chr (sort keys %positions) {
        for my $pos (sort keys %{$positions{$chr}}) {
            my $start = $pos -1;
            $sites_output_fh->print("$chr\t$pos\n");
        }
    }
    $sites_output_fh->close();
    return $segregating_sites_filename;
}

#FIXME Pretty hacky... be less hacky
sub polymutt_dir_for_build {
    my ($self, $build) = @_;
    my $glob_string = $build->data_directory . "/variants/snv/polymutt-*";

    my @all_dirs = glob($glob_string);
    # Ignore -broken links
    my @dirs = grep {! m/broken$/} @all_dirs;

    unless (scalar(@dirs) == 1) {
        die $self->error_message("Found " . scalar(@dirs) . " possible polymutt dirs in " . $build->data_directory . " and I expected to find one.");
    }
    my $dir = $dirs[0];

    unless (-d $dir) {
        die $self->error_message("Polymutt dir $dir does not exist or is not a directory");
    }
    return $dir;
}

# Return the variant vcf (from the model group) for this build.
# This will be the original vcf if region limiting was not requested. If it was requested, return the post-region-limiting file
sub variant_vcf_for_build {
    my ($self, $build, $verify_existance) = @_;
    if ($self->roi_file) {
        return $self->region_limited_variant_vcf_for_build($build, $verify_existance);
    } else {
        return $self->original_vcf_for_build($build);
    }
}

sub region_limited_variant_vcf_for_build {
    my ($self, $build, $verify_existance) = @_;
    my $vcf = $self->subdir_for_build($build) . "/snvs.region_limited.vcf.gz"; 
    if ($verify_existance && !(-s $vcf)) {
        die $self->error_message("Region limited variant vcf $vcf does not exist or has no size");
    }
    return $vcf;
}

sub original_vcf_for_build {
    my ($self, $build) = @_;
    my $vcf = $build->data_directory . "/variants/snvs.vcf.gz";
    unless (-s $vcf) {
        die $self->error_message("Original vcf $vcf does not exist or has no size");
    }
    return $vcf;
}

sub new_glf_index_for_build {
    my ($self, $build) = @_;
    return $self->subdir_for_build($build) . "/polymutt.glfindex"
}

sub old_glf_index_for_build {
    my ($self, $build) = @_;
    my $result = join('/', $self->polymutt_dir_for_build($build), 'polymutt.glfindex');
    unless(-s $result) {
        die $self->error_message("Original glf file $result does not exist or has no size.");
    }
    return $result;
}

sub dat_file_for_build {
    my ($self, $build) = @_;
    my $file = $self->polymutt_dir_for_build($build) . "/polymutt.dat";
    unless (-s $file) {
        die $self->error_message("Dat file $file does not exist or has no size");
    }
    return $file;
}

# Return all glfs contained in the polymutt dir of this build of PhenotypeCorrelation
sub glfs_for_build {
    my ($self, $build) = @_;
    my $polymutt_dir = $self->polymutt_dir_for_build($build);
    my @input_glfs = glob("$polymutt_dir/*.glf");
    unless (@input_glfs) {
        die $self->error_message("No glfs found in $polymutt_dir");
    }
    return @input_glfs;
}

# Return the family id contained in this build of PhenotypeCorrelation
sub family_id_for_build {
    my ($self, $build) = @_;
    my $ped_file = basename($build->pedigree_file_path->path);
    my ($family_id) = $ped_file =~ m/(\w+)\.ped/;
    unless ($family_id) {
        die $self->error_message("Failed to parse out a family id from ped file $ped_file");
    }
    return $family_id;
}

# Return the subdirectory for (new) work done with this build of PhenotypeCorrelation
sub subdir_for_build {
    my ($self, $build) = @_;
    my $subdir = $self->output_dir . "/" . $self->family_id_for_build($build);
    return ($subdir);
}

# Return the force-genotype vcf that has been created for this build
sub force_genotype_vcf_for_build {
    my ($self, $build, $verify_existance) = @_;
    my $vcf = $self->subdir_for_build($build) . "/snvs.vcf.gz";
    if ($verify_existance && !(-s $vcf)) {
        die $self->error_message("Force genotyped vcf $vcf does not exist or has no size");
    }
    return $vcf;
}

# Return the backfilled vcf which was created from the force genotype vcf and the original vcf combined
sub backfilled_vcf_for_build {
    my ($self, $build, $verify_existance) = @_;
    my $vcf = $self->subdir_for_build($build) . "/snvs.backfilled.vcf.gz";
    if ($verify_existance && !(-s $vcf)) {
        die $self->error_message("Backfilled vcf $vcf does not exist or has no size");
    }
    return $vcf;
}

# Run polymutt (just in standard mode, no denovo) with the segregating sites file to force genotype
sub run_polymutt_on_segregating_sites {
    my ($self, $segregating_sites_file) = @_;

    # Params for every operation
    my %params = (
        version => $self->polymutt_version,
        denovo => 0,
        roi_file => $segregating_sites_file,
        bgzip => 1,
        chr2process => $self->chr2process,
        fix_alt_and_gt => 1,
    );

    my @global_input_properties = keys %params;

    my $workflow_model = Workflow::Model->create(
        name => 'Polymutt for segregating sites',
        input_properties => [
            @global_input_properties,
        ],
    );
    $workflow_model->log_dir($self->output_dir);

    my @builds = $self->_builds;
    for my $build (@builds) {
        my $output_vcf = $self->force_genotype_vcf_for_build($build);
        my $ped_file = $build->pedigree_file_path->path;
        my $family_id = $self->family_id_for_build($build);

        # params that differ per operation
        my %params_for_family;
        $params_for_family{$family_id."_output_vcf"} = $output_vcf;
        $params_for_family{$family_id."_glf_index"} = $self->new_glf_index_for_build($build);
        $params_for_family{$family_id."_dat_file"} = $self->dat_file_for_build($build);
        $params_for_family{$family_id."_ped_file"} = $ped_file;
        %params = (%params, %params_for_family);

        # Add the input properties
        my @new_input_connector_properties = keys %params_for_family;
        my $input_connector = $workflow_model->get_input_connector;
        my $input_connector_properties = $input_connector->operation_type->output_properties;
        push @{$input_connector_properties}, @new_input_connector_properties;
        $input_connector->operation_type->output_properties($input_connector_properties);

        my $polymutt_operation = $workflow_model->add_operation(
            name => join(" ",("polymutt", $family_id)),
            operation_type => Workflow::OperationType::Command->get("Genome::Model::Tools::Relationship::RunPolymutt"),
        );

        # Connect global properties
        for my $property (@global_input_properties) {
            $workflow_model->add_link(
                left_operation => $input_connector,
                left_property => $property,
                right_operation => $polymutt_operation,
                right_property => $property,
            );
        }

        # Connect per-family properties
        for my $property (@new_input_connector_properties) {
            my ($right_property) = $property =~ m/$family_id\_(.*)/;
            $workflow_model->add_link(
                left_operation => $input_connector,
                left_property => $property,
                right_operation => $polymutt_operation,
                right_property => $right_property,
            );
        }

        my $new_output_connector_property = "$family_id" . "_output";
        # Add the output properties
        my $output_connector = $workflow_model->get_output_connector;
        my $output_connector_properties = $output_connector->operation_type->input_properties;
        push @{$output_connector_properties}, $new_output_connector_property;
        $input_connector->operation_type->input_properties($output_connector_properties);

        $workflow_model->add_link(
            left_operation => $polymutt_operation,
            left_property => "output_vcf",
            right_operation => $output_connector,
            right_property => $new_output_connector_property,
        );
    }

    # Validate and run the workflow
    my @errors = $workflow_model->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }

    # Dump the workflow
    my $xml_file = Genome::Sys->open_file_for_writing($self->output_dir . "/run_polymutt_on_segregating_sites.xml");
    $workflow_model->save_to_xml(OutputFile => $xml_file);
    $xml_file->close();

    my $result = Workflow::Simple::run_workflow_lsf($workflow_model, %params);
    unless($result){
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("Workflow did not return correctly.");
    }

    return 1;
}

# Combine each force-genotyped vcf with its original vcf to create one "backfilled" vcf for each original build
sub combine_individual_vcfs {
    my $self = shift;

    # Create the necessary strategy file for these merge operations. Maintain original variant call INFO field values.
    my $merge_strategy_file = $self->output_dir . "/joinx_merge_strategy.individual";
    my $strategy_fh = Genome::Sys->open_file_for_writing($merge_strategy_file);
    $strategy_fh->print("default=first\n");
    $strategy_fh->close;

    # Params for every operation
    my %params = (
        use_version => $self->joinx_version,
        sample_priority => 'order',
        merge_samples => 1,
        use_bgzip => 1,
        merge_strategy_file => $merge_strategy_file,
    );

    my @global_input_properties = keys %params;
    my $workflow_model = Workflow::Model->create(
        name => 'Backfill VCF',
        input_properties => [
        @global_input_properties,
        ],
    );
    $workflow_model->log_dir($self->output_dir);

    my @builds = $self->_builds;
    for my $build (@builds) {
        my $force_genotype_vcf = $self->force_genotype_vcf_for_build($build, 1);
        my $original_vcf = $self->variant_vcf_for_build($build);
        my $output_file = $self->backfilled_vcf_for_build($build);
        my $family_id = $self->family_id_for_build($build);

        # params that differ per operation
        # Input files, output file
        my %params_for_family;
        $params_for_family{$family_id."_input_files"} = [$original_vcf, $force_genotype_vcf];
        $params_for_family{$family_id."_output_file"} = $output_file;
        %params = (%params, %params_for_family);

        # Add the input properties
        my @new_input_connector_properties = keys %params_for_family;
        my $input_connector = $workflow_model->get_input_connector;
        my $input_connector_properties = $input_connector->operation_type->output_properties;
        push @{$input_connector_properties}, @new_input_connector_properties;
        $input_connector->operation_type->output_properties($input_connector_properties);

        my $joinx_operation = $workflow_model->add_operation(
            name => join(" ",("joinx-backfill", $family_id)),
            operation_type => Workflow::OperationType::Command->get("Genome::Model::Tools::Joinx::VcfMerge"),
        );

        # Connect global properties
        for my $property (@global_input_properties) {
            $workflow_model->add_link(
                left_operation => $input_connector,
                left_property => $property,
                right_operation => $joinx_operation,
                right_property => $property,
            );
        }

        # Connect per-family properties
        for my $property (@new_input_connector_properties) {
            my ($right_property) = $property =~ m/$family_id\_(.*)/;
            $workflow_model->add_link(
                left_operation => $input_connector,
                left_property => $property,
                right_operation => $joinx_operation,
                right_property => $right_property,
            );
        }

        my $new_output_connector_property = "$family_id" . "_output";
        # Add the output properties
        my $output_connector = $workflow_model->get_output_connector;
        my $output_connector_properties = $output_connector->operation_type->input_properties;
        push @{$output_connector_properties}, $new_output_connector_property;
        $input_connector->operation_type->input_properties($output_connector_properties);

        $workflow_model->add_link(
            left_operation => $joinx_operation,
            left_property => "output_file",
            right_operation => $output_connector,
            right_property => $new_output_connector_property,
        );
    }

    # Validate and run the workflow
    my @errors = $workflow_model->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }

    # Dump the workflow
    my $xml_file = Genome::Sys->open_file_for_writing($self->output_dir . "/combine_individual_vcfs.xml");
    $workflow_model->save_to_xml(OutputFile => $xml_file);
    $xml_file->close();

    my $result = Workflow::Simple::run_workflow_lsf($workflow_model, %params);
    unless($result){
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("Workflow did not return correctly.");
    }

    return 1;
}

# Combine all backfilled vcfs into one final vcf
sub create_final_vcf {
    my $self = shift;

    my $output_file = $self->output_dir . "/snvs.vcf.gz";
    my @builds = $self->_builds;
    my @backfilled_vcfs;
    for my $build (@builds) {
        my $per_build_output_file = $self->backfilled_vcf_for_build($build, 1);
        push @backfilled_vcfs, $per_build_output_file;
    }

    my $working_directory = $self->output_dir . "/sub_merge_working_dir";
    Genome::Sys->create_directory($working_directory);

    # Create the necessary strategy file for this merge operation. Ignore all INFO fields since joining them makes little sense cross sample.
    my $merge_strategy_file = $self->output_dir . "/joinx_merge_strategy.final";
    my $strategy_fh = Genome::Sys->open_file_for_writing($merge_strategy_file);
    $strategy_fh->print("default=ignore\n");
    $strategy_fh->close;

    my $command = Genome::Model::Tools::Joinx::SafeVcfMerge->create(
        use_version => $self->joinx_version,
        input_files => \@backfilled_vcfs,
        use_bgzip => 1,
        output_file => $output_file,
        max_files_per_merge => 100, # FIXME the default is probably fine, let's test it
        remove_intermediate_files => 0, # FIXME probably stop doing this once debugging is done
        merge_strategy_file => $merge_strategy_file,
        working_directory => $working_directory,
    );

    $self->status_message("Running joinx vcf-merge to join all backfilled vcfs to output $output_file");
    unless ($command->execute) {
        die $self->error_message("Failed to run joinx vcf-merge on all backfilled vcfs to output $output_file");
    }

    return 1;
}

# FIXME just copypasted, fix this up
# Region limit the output files
sub region_limit_model_group {
    my $self = shift;

    unless ($self->roi_file) {
        $self->warning_message("No roi_file set, skipping region limiting. Is this intentional? With no roi set, this may take a LONG time.");
        return 1;
    }

    my @builds = $self->_builds;

    my @answers;
    my %in_out;

    my %inputs;
    $inputs{joinx_version} = $self->joinx_version;

    $inputs{region_bed_file} = $self->roi_file;
    $inputs{roi_name} = $self->roi_name;
    $inputs{wingspan} = $self->wingspan;

    my @inputs;
    my $count=1;

    #set up individualized input params and input values
    for my $build (@builds){
        my $sample = $build->model->subject->name;
        my $vcf = $self->original_vcf_for_build($build);
        my $output = $self->region_limited_variant_vcf_for_build($build);
        $in_out{$vcf} = $output;
        push @inputs, ("input_vcf_".$count,"output_vcf_".$count);
        $inputs{"input_vcf_".$count} = $vcf;
        $inputs{"output_vcf_".$count} = $output;
        push @answers, $output;
        $count++;
    }

    my $workflow = Workflow::Model->create(
        name => 'Multi-Vcf Merge',
        input_properties => [
            "region_bed_file",
            "roi_name",
            "wingspan",
            "joinx_version",
            @inputs,
        ],
        output_properties => [
            'output',
        ],
    );

    $workflow->log_dir($self->output_dir);

    #add individual region-limiting operations
    for my $num (1..($count-1)){
        my $region_limit_operation = $workflow->add_operation(
            name => "region limiting ".$num,
            operation_type => Workflow::OperationType::Command->get("Genome::Model::Tools::Vcf::RegionLimit"),
        );

        #link common properties
        for my $prop ("region_bed_file","wingspan","roi_name"){
            $workflow->add_link(
                left_operation => $workflow->get_input_connector,
                left_property => $prop,
                right_operation => $region_limit_operation,
                right_property => $prop,
            );
        }

        #link individual inputs and outputs
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => "input_vcf_".$num,
            right_operation => $region_limit_operation,
            right_property => "vcf_file",
        );
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => "output_vcf_".$num,
            right_operation => $region_limit_operation,
            right_property => "output_file",
        );

        #link to output
        $workflow->add_link(
            left_operation => $region_limit_operation,
            left_property => "output_file",
            right_operation => $workflow->get_output_connector,
            right_property => "output",
        );
    }

    #validate workflow
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating region-limiting workflow\n";
    }

    # Dump the workflow
    my $xml_file = Genome::Sys->open_file_for_writing($self->output_dir . "/region_limiting_workflow.xml");
    $workflow->save_to_xml(OutputFile => $xml_file);
    $xml_file->close();

    $self->status_message("Now launching the region-limiting workflow.");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);

    unless($result){
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("Workflow did not return correctly.");
    }

    #check output files to make sure they exist
    if(my @error = grep{ not(-e $_)} @answers){
        die $self->error_message("The following region limit output files could not be found: ".join("\n",@error));
    }

    #return a list of the output files
    return @answers;
}

1;
