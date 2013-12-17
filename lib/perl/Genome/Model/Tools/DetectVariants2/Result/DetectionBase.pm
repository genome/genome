package Genome::Model::Tools::DetectVariants2::Result::DetectionBase;

use strict;
use warnings;

use Sys::Hostname;
use File::Path 'rmtree';

use Genome;

class Genome::Model::Tools::DetectVariants2::Result::DetectionBase {
    is => ['Genome::Model::Tools::DetectVariants2::Result::Base'],
    is_abstract => 1,
    sub_classification_method_name => '_resolve_subclass_name',
    has => [
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'reference_build_id',
        },
        region_of_interest => {
            is => 'Genome::FeatureList',
            is_optional => 1,
            id_by => 'region_of_interest_id',
        },
        _disk_allocation => {
            is => 'Genome::Disk::Allocation',
            is_optional => 1,
            is_many => 1,
            reverse_as => 'owner'
        },
    ],
    has_param => [
        detector_name => {
            is => 'Text',
            doc => 'The name of the detector to use',
        },
        detector_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'Additional parameters to pass to the detector',
        },
        detector_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'Version of the detector to use',
        },
        chromosome_list => {
            is => 'Text',
            is_optional => 1,
            doc => 'The chromosome(s) on which the detection was run',
        },
    ],
    has_optional_input => [
        # old API
        aligned_reads => {
            is => 'Text',
            is_optional => 1,
            doc => 'The path to the aligned reads file',
        },
        control_aligned_reads => {
            is => 'Text',
            is_optional => 1,
            doc => 'The path to the control aligned reads file',
        },
        # new API
        alignment_results => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged',
            is_optional => 1,
            is_many => 1,
            doc => 'The path to the aligned reads file',
        },
        control_alignment_results => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged',
            is_optional => 1,
            is_many => 1,
            doc => 'The path to the control aligned reads file',
        },
        roi_list => {
            is => 'Genome::FeatureList',
            is_optional => 1,
            doc => 'only variants in these regions will be included in the final VCF',
        },
        roi_wingspan => {
            is => 'Number',
            is_optional => 1,
            doc => 'include variants within N nucleotides of a region of interest'
        },
        pedigree_file_path => {
            is => 'FilePath',
            is_optional => 1,
            doc => 'when supplied overrides the automatic lookup of familial relationships'
        },
    ],
    has_input => [
        reference_build_id => {
            is => 'Number',
            doc => 'the reference to use by id',
        },
        region_of_interest_id => {
            is => 'Text',
            doc => 'The feature-list representing the region of interest (if present, only variants in the set will be reported)',
            is_optional => 1,
        },
    ],
    has_transient_optional => [
        _instance => {
            is => 'Genome::Model::Tools::DetectVariants2::Base',
            doc => 'The instance of the entity that is creating this result',
        }
    ],
    doc => 'This class represents the result of a detect-variants operation.',
};

sub create {
    my $class = shift;

    #This will do some locking and the like for us.
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    eval {
        $self->_prepare_output_directory;

        my $instance = $self->_instance;
        my $instance_output = $instance->output_directory;
        $self->_check_instance_output($instance_output);

        Genome::Sys->create_symlink_and_log_change($self, $self->output_dir, $instance_output);

        $instance->_generate_result;
        $self->_set_result_file_permissions;
    };
    if($@) {
        my $error = $@;
        $self->_cleanup;
        die $error;
    }

    $self->status_message("Resizing the disk allocation...");
    if ($self->_disk_allocation) {
        my $result = eval { $self->_disk_allocation->reallocate };
        if($@ or not $result) {
            my $err = $@;
            $self->warning_message("Failed to reallocate disk allocation: " . $self->_disk_allocation->id . ($err || ''));
        }
    }

    return $self;
}

sub _check_instance_output {
    my ($class, $instance_output) = @_;
    $class = ref $class if ref $class;

    if (-l $instance_output or -e $instance_output) {
        # If the detector output is not a symlink, it was generated before these were software results.
        # Archive the existing stuff and regenerate so we get a nifty software result.
        if (not -l $instance_output and -d $instance_output) {
            my ($parent_dir, $sub_dir) = $instance_output =~ /(.+)\/(.+)/;
            die $class->error_message("Unable to determine parent directory from $instance_output.") unless $parent_dir;
            die $class->error_message("Unable to determine sub-directory from $instance_output.") unless $sub_dir;
            die $class->error_message("Parse error when determining parent directory and sub-directory from $instance_output.") unless ($instance_output eq "$parent_dir/$sub_dir");

            my $archive_name = "old_$sub_dir.tar.gz";
            die $class->error_message("Archive already exists, $parent_dir/$archive_name.") if (-e "$parent_dir/$archive_name");

            $class->status_message('Archiving old non-software-result ' . $instance_output . " to $archive_name.");
            system("cd $parent_dir && tar -zcvf $archive_name $sub_dir && rm -rf $sub_dir");
        }
        # In this case, a symlink exists...
        else {
            $class->warning_message('Instance output directory (' . $instance_output . ') already exists!');
            my $allocation_dir = readlink $instance_output;
            my @parts = split "-", $allocation_dir;
            my $allocation_owner_id = $parts[-1];
            my $result = Genome::SoftwareResult->get($allocation_owner_id);
            my $allocation = Genome::Disk::Allocation->get(owner_id => $allocation_owner_id);

            if ($allocation and not $result) {
                # Allocation exists without a result the whole time the result is being created. Ideally locks
                # would prevent us from getting here during that window but our locks are not 100% reliable.
                my @error_message = (
                    "Found allocation at ($allocation_dir) but no software result for it's owner ID ($allocation_owner_id).",
                    "This is either because the software result is currently being generated or because the allocation has been orphaned.",
                    "If it is determined that the allocation has been orphaned then the allocation will need to be removed.",
                );
                die $class->error_message(join(' ', @error_message));
            }
            # If a test name is set, we can remove the symlink and proceed
            elsif ($result and defined $result->test_name and -l $instance_output) {
                $class->warning_message("The software result for the existing symlink has a test name set; removing symlink.");
                unlink $instance_output;
            }
            elsif ($result and not $allocation) {
                # A result without an allocation... this really shouldn't ever happen, unless someone deleted the allocation row from the database?
                die $class->error_message("Found a software result (" . $result->__display_name__ . ") that has output directory " .
                    "($instance_output) but no allocation.");
            }
            elsif ($result and $allocation) {
                # Finding a result and an allocation means either:
                # 1) This work was already done, but for whatever reason we didn't find the software result before we decided to do the work.
                # 2) We're doing different work but pointing it at a place where work has already been done for something else. Can't replace it.
                die $class->error_message("Found allocation and software result for path $instance_output, cannot create new result!");
            }
            elsif (!$result && !$allocation && -l $instance_output && ! -e $allocation_dir) {
                $class->warning_message("No allocation or software result and symlink ($instance_output) target ($allocation_dir) does not exist; removing symlink.");
                unlink $instance_output;
            }
            else {
                die $class->error_message("Unexpected condition in " . __PACKAGE__ . "::_check_instance_output");
            }
        }
    }

    return 1;
}

sub _gather_params_for_get_or_create {
    my $class = shift;

    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }
    }

    #my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    #my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my $subclass_name = (exists $is_param{filter_name} ? 'Genome::Model::Tools::DetectVariants2::Result::Filter' : 'Genome::Model::Tools::DetectVariants2::Result');
    if(exists($is_param{vcf_version})){
        $subclass_name = 'Genome::Model::Tools::DetectVariants2::Result::DetectorVcf';
    }

    my %software_result_params = (#software_version=>$params_bx->value_for('aligner_version'),
        #params_id=>$params_bx->id,
        #inputs_id=>$inputs_bx->id,
                                  subclass_name=>$subclass_name);

    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs=>\%is_input,
        params=>\%is_param,
        _instance => (exists $params{instance} ? $params{instance} : undef),
    };
}

sub estimated_kb_usage {
    my $self = shift;

    return 10_000_000; #TODO be more dynamic about this
}

#TODO standardize this with _resolve_allocation_subdirectory in other DV2 Results and pull up to Result::Base
sub _resolve_subdirectory {
    my $self = shift;

    my $hostname = hostname;
    my $user = $ENV{'USER'};
    my $base_dir = sprintf("detect-variants--%s-%s-%s-%s", $hostname, $user, $$, $self->id);
    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments', $base_dir);
    return $directory;
}

sub _prepare_output_directory {
    my $self = shift;

    return $self->output_dir if $self->output_dir;

    my $subdir = $self->_resolve_subdirectory;
    unless ($subdir) {
        $self->error_message("failed to resolve subdirectory for instrument data.  cannot proceed.");
        die $self->error_message;
    }

    my $allocation = $self->_disk_allocation;

    unless($allocation) {
        my %allocation_parameters = (
            disk_group_name => $ENV{GENOME_DISK_GROUP_MODELS},
            allocation_path => $subdir,
            owner_class_name => $self->class,
            owner_id => $self->id,
            kilobytes_requested => $self->estimated_kb_usage,
        );

        $allocation = Genome::Disk::Allocation->allocate(%allocation_parameters);
    }

    my $output_dir = $allocation->absolute_path;
    unless (-d $output_dir) {
        $self->error_message("Allocation path $output_dir doesn't exist!");
        die $self->error_message;
    }

    $self->output_dir($output_dir);

    return $output_dir;
}

sub _cleanup {
    my $self = shift;

    my $instance = $self->_instance;
    if($instance) {
        my $instance_output = $instance->output_directory;
        # Remove trailing slashes so readlink will work
        if ($instance_output =~ m/\/+$/) {
            $instance_output =~ s/\/+$//;
        }
        if(readlink($instance_output) eq $self->output_dir) {
            unlink($instance_output);
        }
    }

    return unless $self->_disk_allocation;

    $self->status_message('Now deleting allocation with owner_id = ' . $self->id);
    my $allocation = $self->_disk_allocation;
    $allocation->deallocate if $allocation;
}

sub _resolve_subclass_name {
    my $class = shift;

    if (ref($_[0]) and $_[0]->isa(__PACKAGE__)) {
        my $filter_name = $_[0]->params(name => 'filter_name');
        return $filter_name ? 'Genome::Model::Tools::DetectVariants2::Result::Filter' : 'Genome::Model::Tools::DetectVariants2::Result';
    }
    else {
        my $filter_name = $class->define_boolexpr(@_)->value_for('filter_name');
        return $filter_name ? 'Genome::Model::Tools::DetectVariants2::Result::Filter' : 'Genome::Model::Tools::DetectVariants2::Result';
    }
    return;
}

sub _set_result_file_permissions {
    my $self = shift;
    my $output_dir = $self->output_dir;
    if($output_dir =~ m/\/$/){
        $output_dir =~ s/\/$//;
    }

    chmod 02775, $output_dir;
    for my $subdir (grep { -d $_  } glob("$output_dir/*")) {
        chmod 02775, $subdir;
    }

    # Make everything in here read-only
    for my $file (grep { -f $_  } glob("$output_dir/*")) {
        chmod 0444, $file;
    }
}


1;
