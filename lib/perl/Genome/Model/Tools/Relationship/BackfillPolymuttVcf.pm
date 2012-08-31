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
    ],
    has_optional_input => [
       roi_file => {
           is => "Text",
           doc => "Optional roi file (just chrom pos) to use instead of generating one from segregating sites in the model group",
       },
       polymutt_version=> {
           is => "Text",
           default => 'vcf.0.01',
       },
       joinx_version=> {
           is => "Text",
           default => '1.6',
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
            default => 'long',
        },
    ],
};

sub help_detail {
    'Use polymutt to force-genotype all segregating sites and produce one final vcf'
}

#/gscuser/dlarson/src/polymutt.0.01/bin/polymutt -p 20000492.ped -d 20000492.dat -g 20000492.glfindex --minMapQuality 1 --nthreads 4 --vcf 20000492.standard.vcf
sub execute {
    my $self=shift;

    $self->status_message("Resolving valid builds for the given model group");
    $self->resolve_valid_builds;

    $self->status_message("Preparing directories");
    $self->prepare_directories;

    $self->status_message("Resolving the roi file to be used");
    my $roi_file = $self->resolve_roi_file;

    $self->status_message("Running Polymutt with the roi file");
    $self->run_polymutt_on_roi($roi_file);

    $self->status_message("Backfilling all original vcfs with force-genotype data");
    $self->combine_individual_vcfs;

    $self->status_message("Combining all backfilled vcfs into one final vcf");
    $self->create_final_vcf;

    return 1;
}

sub resolve_valid_builds {
    my $self = shift;
    my $mg = $self->model_group;
    my @builds;

    for my $model ($mg->models) {
        my $build = $model->current_build;
        if($build && $build->status eq 'Succeeded' && $build->qc_succeeded) {
            push @builds, $build;
        } else {
            $self->status_message("Skipping " . $model->name . "No successful builds that passed QC");
        }
    }

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

sub resolve_roi_file {
    my $self = shift;

    my $roi_file = $self->roi_file;
    unless($roi_file) {
        $self->status_message("No roi file supplied, attempting to make a roi file of all possible sites from Model group outputs...");
        $roi_file = $self->assemble_list_of_segregating_sites();
    }
    unless (-s $roi_file) {
        die $self->error_message("roi file $roi_file does not exist or is empty");
    }

    return $roi_file;
}

# Gather a list of all variant sites from all builds
# FIXME can probably be replaced with a call to joinx vcf-merge piped to cut -f1,2 to produce chr pos
sub assemble_list_of_segregating_sites {
    my $self = shift;

    my $output_dir = $self->output_dir;
    my @builds = $self->_builds;
    my $roi_filename  = $output_dir . "/segregating_sites.bed";

    my %positions;
    #open every final original snvs.vcf.gz and accumulate all unique positions that will need to be polled across the model group. (CHR POS)
    for my $build (@builds) {
        my $snvs_vcf = $self->original_vcf_for_build($build);
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

    #print roi file
    my $roi_output_fh = Genome::Sys->open_file_for_writing($roi_filename);
    for my $chr (sort keys %positions) {
        for my $pos (sort keys %{$positions{$chr}}) {
            my $start = $pos -1;
            $roi_output_fh->print("$chr\t$pos\n");
        }
    }
    $roi_output_fh->close();
    return $roi_filename;
}

#FIXME Pretty hacky... be less hacky
sub polymutt_dir_for_build {
    my ($self, $build) = @_;
    my $dir = $build->data_directory . "/variants/snv/polymutt-0.11-8b1c96ee44b70a5c936f48fd06d74d07";
    unless (-d $dir) {
        die $self->error_message("Polymutt dir $dir does not exist or is not a directory");
    }
    return $dir;
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

# Return the force-genotype vcf (before header fixing) that has been created for this build
sub pre_fix_force_genotype_vcf_for_build {
    my ($self, $build, $verify_existance) = @_;
    my $vcf = $self->subdir_for_build($build) . "/pre_vcf_fix.snvs.vcf.gz";
    if ($verify_existance && !(-s $vcf)) {
        die $self->error_message("Pre-vcf-fix force genotyped vcf $vcf does not exist or has no size");
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

# Run polymutt (just in standard mode, no denovo) with the roi file to force genotype
sub run_polymutt_on_roi {
    my ($self, $roi_file) = @_;

    # Params for every operation
    my %params = (
        version => $self->polymutt_version,
        denovo => 0,
        roi_file => $roi_file,
        bgzip => 1,
        chr2process => $self->chr2process,
    );

    my @global_input_properties = keys %params;

    my $workflow_model = Workflow::Model->create(
        name => 'Polymutt for ROI',
        input_properties => [
            @global_input_properties,
        ],
    );
    $workflow_model->log_dir($self->output_dir);

    my @builds = $self->_builds;
    for my $build (@builds) {
        my $pre_fix_vcf = $self->pre_fix_force_genotype_vcf_for_build($build);
        my $ped_file = $build->pedigree_file_path->path;
        my $family_id = $self->family_id_for_build($build);

        # params that differ per operation
        my %params_for_family;
        $params_for_family{$family_id."_output_vcf"} = $pre_fix_vcf;
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
    my $result = Workflow::Simple::run_workflow_lsf($workflow_model, %params);
    unless($result){
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("Workflow did not return correctly.");
    }

    # Fix up the vcfs
    # FIXME move this to RunPolymutt, or the workflow, rather than doing it inline
    for my $build (@builds) {
        my $pre_fix_vcf = $self->pre_fix_force_genotype_vcf_for_build($build);
        my $post_fix_vcf = $self->force_genotype_vcf_for_build($build);
        $self->status_message("Fixing polymutt vcf $pre_fix_vcf and storing the fix at $post_fix_vcf");
        $self->write_fixed_vcf($pre_fix_vcf, $post_fix_vcf);
    }

    return 1;

}


# Combine each force-genotyped vcf with its original vcf to create one "backfilled" vcf for each original build
sub combine_individual_vcfs {
    my $self = shift;

    # Params for every operation
    my %params = (
        use_version => $self->joinx_version,
        sample_priority => 'order',
        merge_samples => 1,
        use_bgzip => 1,
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
        my $original_vcf = $self->original_vcf_for_build($build);
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
        my $per_build_output_file = $self->backfilled_vcf_for_build($build);
        push @backfilled_vcfs, $per_build_output_file;
    }

    my $command = Genome::Model::Tools::Joinx::VcfMerge->create(
        use_version => $self->joinx_version,
        input_files => \@backfilled_vcfs,
        use_bgzip => 1,
        output_file => $output_file,
    );

    $self->status_message("Running joinx vcf-merge to join all backfilled vcfs to output $output_file");
    unless ($command->execute) {
        die $self->error_message("Failed to run joinx vcf-merge on all backfilled vcfs to output $output_file");
    }

    return 1;
}

# FIXME copied... call this in Genome::Model::Tools::Relationship::MergeAndFixVcfs or move to a base class
# Copy is a conglomoration of write_merged_header and fix_vcf_header
sub write_fixed_vcf {
    my ($self, $input_vcf, $output_vcf) = @_;
    my @info_lines_to_add = (qq|##INFO=<ID=DQ,Number=1,Type=Float,Description="De Novo Mutation Quality">|);
    push @info_lines_to_add, qq|##INFO=<ID=DA,Number=1,Type=Integer,Description="De Novo Mutation Allele">|;
    my @format_lines_to_add = (qq|##FORMAT=<ID=DNGL,Number=10,Type=Integer,Description="Denovo Genotype Likelihoods">|);
    push @format_lines_to_add, qq|##FORMAT=<ID=DNGT,Number=1,Type=String,Description="Genotype">|;
    push @format_lines_to_add, qq|##FORMAT=<ID=DNGQ,Number=1,Type=Integer,Description="Genotype Quality">|;
    ####hack to add in header for bingshan's tag that he hasn't added himself
    my @missing_info_lines = $self->missing_info_lines($input_vcf);
    push @info_lines_to_add, @missing_info_lines if (@missing_info_lines);
    ######

    my $fh = Genome::Sys->open_gzip_file_for_writing($output_vcf);
    my $ifh = Genome::Sys->open_gzip_file_for_reading($input_vcf);
    my ($info_printed, $format_printed)=(0,0);
    while(my $line = $ifh->getline) {
        if($line =~m/^#/) {
            # Fix incorrect data types
            if($line =~m/ID=PS/) {
                $line =~ s/Integer/Float/;
            }
            #Fix spacing before description
            $line =~s/, Description/,Description/;

            # Fix invalid type labels
            $line =~s/,String/,Type=String/;

            # Fix lines without a "number" tag
            if ( ($line =~ m/^##INFO/ or $line =~ m/^##FORMAT/) and not ($line =~ m/Number=/) ) {
                $line =~ s/,/,Number=1,/;
            }

            if($line =~m/#INFO/ && !$info_printed) {
                for my $info_line (@info_lines_to_add) {
                    $fh->print($info_line ."\n");
                }
                $info_printed=1;
            }
            if($line =~m/#FORMAT/ && !$format_printed) {
                for my $format_line (@format_lines_to_add) {
                    $fh->print($format_line ."\n");
                }
                $format_printed=1
            }
            $fh->print($line);
        }
        else {
            $fh->print($line);
        }
    }
    $fh->close; $ifh->close;

    return 1;
}

# FIXME copied... call this in Genome::Model::Tools::Relationship::MergeAndFixVcfs or move to a base class
sub missing_info_lines {
    my ($self, $input_vcf) = @_;

    my %possible_tags = (
        AB => qq|##INFO=<ID=AB,Number=1,Type=Float,Description="Allelic Balance">|,
        BA => qq|##INFO=<ID=BA,Number=1,Type=String,Description="Best Alternative Allele">|,
    );

    my @info_lines_to_add;
    for my $tag (keys %possible_tags) {
        $DB::single=1;
        chomp(my @number_of_tags = `zcat $input_vcf | grep $tag`);
        if(@number_of_tags) {
            my $has_header=0;
            for my $line (@number_of_tags) {
                if($line=~m/^##INFO=<ID=$tag,/) {
                    $has_header=1;
                }
            }
            unless($has_header) {
                push @info_lines_to_add, $possible_tags{$tag};
            }
        }
    }

    return @info_lines_to_add;
}

1;
