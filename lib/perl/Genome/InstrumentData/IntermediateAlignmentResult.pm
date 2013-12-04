package Genome::InstrumentData::IntermediateAlignmentResult;

use Data::Dumper;
use Genome;
use Sys::Hostname;
use File::Basename;
use Carp qw/confess/;

use warnings;
use strict;

class Genome::InstrumentData::IntermediateAlignmentResult {
    is_abstract => 1,
    is=>['Genome::SoftwareResult::Stageable'],
    sub_classification_method_name => '_resolve_subclass_name',   
    has => [
        instrument_data         => {
                                    is => 'Genome::InstrumentData',
                                    id_by => 'instrument_data_id'
                                },
        aligner_index           => {
                                    is => 'Genome::Model::Build::ReferenceSequence::AlignerIndex',
                                    id_by => 'aligner_index_id',
                                },
        reference_name          => { via => 'aligner_index', to => 'reference_name', is_mutable => 0, is_optional => 1 },

        aligner                 => { 
                                    calculate_from => [qw/aligner_name aligner_version aligner_params/], 
                                    calculate => q|no warnings; "$aligner_name $aligner_version $aligner_params"| 
                                },
        
        _disk_allocation        => { is => 'Genome::Disk::Allocation', is_optional => 1, is_many => 1, reverse_as => 'owner' },

    ],
    has_input => [
        input_file              => {
                                    is=>'Text',
                                    doc=>'Path to the (potentially pre-processed) input fastq or bam file',
                                },
        instrument_data_id      => {
                                    is => 'Number',
                                    doc => 'the local database id of the instrument data (reads) to align',
                                },
        instrument_data_segment_type => {
                                    is => 'String',
                                    doc => 'Type of instrument data segment to limit within the instrument data being aligned (e.g. "read_group")',
                                    is_optional => 1,
        },
        instrument_data_segment_id => {
                                    is => 'String',
                                    doc => 'Identifier for instrument data segment to limit within the instrument data being aligned (e.g. read group ID)',
                                    is_optional => 1,
        },
        aligner_index_id        => {
                                    is => 'Number',
                                    doc => 'the aligner index to use by id',
                                },
    ],
    has_param => [
        input_pass              => {
                                    is => 'Integer',
                                    doc => 'Specify 1 or 2 to indicate which set of reads to use for paired end reads in .bam files.',
                                    is_optional => 1,
                                },
        aligner_name            => {
                                    is => 'Text', default_value => 'maq',
                                    doc => 'the name of the aligner to use, maq, blat, newbler etc.',
                                },
        aligner_version         => {
                                    is => 'Text',
                                    doc => 'the version of the aligner to use, i.e. 0.6.8, 0.7.1, etc.',
                                    is_optional=>1,
                                },
        aligner_params          => {
                                    is => 'Text',
                                    is_optional=>1,
                                    doc => 'any additional params for the aligner',
                                },
        force_fragment          => {
                                    is => 'Boolean',    
                                    is_optional=>1,
                                    doc => 'Force this run to be treated as a fragment run, do not do pairing',
                                },
        samtools_version        => {
                                    is => 'Text',
                                    is_optional => 1,
                                    doc => 'the samtools version used',
                                },
        trimmer_name            => {
                                    is => 'Text',
                                    doc => 'Trimmer strategy used to create input file.',
                                    is_optional=>1,
                                },
        trimmer_version         => {
                                    is => 'Text',
                                    doc => 'Trimmer version to used to create input file.',
                                    is_optional=>1,
                                },
        trimmer_params          => {
                                    is => 'Text',
                                    is_optional=>1,
                                    doc => 'Trimmer parameters used to create input file.',
                                },
    ],
    has_transient => [
        temp_scratch_directory  => {
                                    is=>'Text',
                                    doc=>'Temp scratch directory',
                                    is_optional=>1,
                                },
        flagstat_file           => {
                                    is =>'Text',
                                    doc => 'Precalculated flagstat file for the BAM (if any)',
                                    is_optional=>1,
                                },
    ],
};

sub _run_aligner {
    confess "unimplemented method: please define _run_aligner in your alignment subclass.";
}

sub required_arch_os {
    # override in subclasses if 64-bit is not required
    'x86_64' 
}

sub required_rusage { 
    # override in subclasses
    # e.x.: "-R 'select[model!=Opteron250 && type==LINUX64] span[hosts=1] rusage[tmp=50000:mem=12000]' -M 1610612736";
    ''
}

sub _working_dir_prefix {
    "alignment";
}

sub _resolve_subclass_name {
    my $class = shift;

    if (ref($_[0]) and $_[0]->isa(__PACKAGE__)) {
        my $aligner_name = $_[0]->aligner_name;
        return join('::', 'Genome::InstrumentData::IntermediateAlignmentResult', $class->_resolve_subclass_name_for_aligner_name($aligner_name));
    }
    elsif (my $aligner_name = $class->define_boolexpr(@_)->value_for('aligner_name')) {
        return join('::', 'Genome::InstrumentData::IntermediateAlignmentResult', $class->_resolve_subclass_name_for_aligner_name($aligner_name));
    }
    return;
}

sub _resolve_subclass_name_for_aligner_name {
    my ($class,$aligner_name) = @_;
    my @type_parts = split(' ',$aligner_name);

    my @sub_parts = map { ucfirst } @type_parts;
    my $subclass = join('',@sub_parts);

    return $subclass;
}

sub verify_os {
    my $class = shift;
    my $actual_os = Genome::Config->arch_os();
    $class->status_message("OS is $actual_os");
    my $required_os = $class->required_arch_os;
    $class->status_message("Required OS is $required_os");
    unless ($required_os eq $actual_os) {
        confess $class->error_message("This logic can only be run on a $required_os machine!  (running on $actual_os)");
    }
}

sub verify_disk_space {
    my $self = shift;

    # TODO: move disk_group, estimated_size, allocation and promotion up into the software result logic
    my $estimated_kb_usage = $self->estimated_kb_usage;
    $self->status_message("Estimated disk for this data set: " . $estimated_kb_usage . " kb");
    $self->status_message("Check for available disk...");
    my @available_volumes = Genome::Disk::Volume->get(disk_group_names => $ENV{GENOME_DISK_GROUP_ALIGNMENTS}); 
    $self->status_message("Found " . scalar(@available_volumes) . " disk volumes");
    my $unallocated_kb = 0;
    for my $volume (@available_volumes) {
        $unallocated_kb += $volume->unallocated_kb;
    }
    $self->status_message("Available disk: " . $unallocated_kb . " kb");
    my $factor = 20;
    unless ($unallocated_kb > ($factor * $estimated_kb_usage)) {
        $self->error_message("NOT ENOUGH DISK SPACE!  This step requires $factor x as much disk as the job will use to be available before starting.");
        confess $self->error_message();
    }
}

sub get {
    my $class = shift;

    # If get is called on a blessed object, we are looking for properties on that instance, not objects.
    # Defining a BoolExpr in that case will cause things to fail.
    return $class->SUPER::get(@_) if ref($class);

    return $class->SUPER::get($class->_process_params(@_));
}

sub _process_params {
    my $class = shift;
    my @params = @_;

    my ($bx, @extra) = $class->define_boolexpr(@params);
    if($bx->specifies_value_for('input_file')) {
        my $absolute_path = $bx->value_for('input_file');
        $bx = $bx->remove_filter('input_file')->add_filter(input_file => basename($absolute_path));
    }

    return ($bx, @extra);
}

sub _modify_params_for_lookup_hash {
    my ($class, $params_ref) = @_;

    my $absolute_path = $params_ref->{'input_file'};
    if (defined($absolute_path)) {
        $params_ref->{'input_file'} = basename($absolute_path);
    }
}

sub create {
    my $class = shift;
    if ($class eq __PACKAGE__ or $class->__meta__->is_abstract) {
        # this class is abstract, and the super-class re-calls the constructor from the correct subclass
        return $class->SUPER::create(@_);
    }

    $class->verify_os();

    #have to process params in advance so UR index is built correctly
    my $self = $class->SUPER::create($class->_process_params(@_));
    return unless $self;

    my ($bx, @extra) = $class->define_boolexpr(@_);
    my $full_input_path = $bx->value_for('input_file');

    if (my $output_dir = $self->output_dir) {
        if (-d $output_dir) {
            $self->status_message("BACKFILL DIRECTORY: $output_dir!");
            return $self;
        }
    }

    $self->verify_disk_space();
    $self->_prepare_working_and_staging_directories;

    # PREPARE THE ALIGNMENT DIRECTORY ON NETWORK DISK
    $self->status_message("Preparing the output directory...");
    $self->status_message("Staging disk usage is " . $self->_staging_disk_usage . " KB");
    my $output_dir = $self->output_dir || $self->_prepare_output_directory;
    $self->status_message("Alignment output path is $output_dir");

    # symlink the input file in the working directory, drop the full path from the input
    my $symlink_target = $self->temp_scratch_directory . "/" . $self->input_file;
    unless(symlink($full_input_path, $symlink_target)) {
        confess "Failed to create symlink of input file " . $self->input_file . " at $symlink_target.";
    }

    # do some work!
    $self->_run_aligner();

    # PROMOTE THE DATA INTO ALIGNMENT DIRECTORY
    $self->status_message("Moving results to network disk...");
    my $product_path;
    unless($product_path = $self->_promote_data) {
        confess $self->error_message("Failed to de-stage data into alignment directory " . $self->error_message);
    }
    
    # RESIZE THE DISK
    # TODO: move this into the actual original allocation so we don't need to do this 
    $self->status_message("Resizing the disk allocation...");
    if ($self->_disk_allocation) {
        my %params;
        $params{allow_reallocate_with_move} = 0;
        $params{allow_reallocate_with_move} = 1 if $self->_disk_allocation->kilobytes_requested < 10_485_760; # 10GB
        unless ($self->_disk_allocation->reallocate(%params)) {
            $self->warning_message("Failed to reallocate my disk allocation: " . $self->_disk_allocation->id);
        }
        $self->output_dir($self->_disk_allocation->absolute_path); #update if was moved
    }
        
    $self->status_message("Intermediate alignment result generation complete.");
    return $self;
}

sub _gather_params_for_get_or_create {
    my $class = shift;
    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my $aligner_name = $bx->value_for('aligner_name');
    my $subclass = join '::', 'Genome::InstrumentData::IntermediateAlignmentResult', $class->_resolve_subclass_name_for_aligner_name($aligner_name);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $subclass->__meta__;
    for my $key ($subclass->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key}; 
        }

    }

    #my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($subclass, %is_input);
    #my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($subclass, %is_param);

    my %software_result_params = (#software_version=>$params_bx->value_for('aligner_version'),
        #params_id=>$params_bx->id,
        #inputs_id=>$inputs_bx->id,
                                  subclass_name=>$subclass);
    
    return {
        software_result_params => \%software_result_params,
        subclass => $subclass,
        inputs=>\%is_input,
        params=>\%is_param,
    };
}

sub _prepare_working_and_staging_directories {
    my $self = shift;

    $self->status_message("Prepare working directories...");

    unless ($self->_prepare_staging_directory) {
        $self->error_message("Failed to prepare staging directory");
        return;
    }
    my $hostname = hostname;
    my $user = $ENV{'USER'};

    my $scratch_tempdir = Genome::Sys->create_temp_directory();
    $self->temp_scratch_directory($scratch_tempdir);
    unless(-d $scratch_tempdir) {
        confess "Failed to create a temp scrach directory for working files.";
    }

    $self->status_message("Staging path is " . $self->temp_staging_directory);
    $self->status_message("Working path is " . $self->temp_scratch_directory);
    return 1;
} 

sub estimated_kb_usage {
    30000000;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $instrument_data = $self->instrument_data;
    my $staged_basename = File::Basename::basename($self->temp_staging_directory);
    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'alignment_data',$instrument_data->id,$staged_basename);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_ALIGNMENTS};
}

# this behavior was in the alignment class earlier and was set as changeable in SoftwareResult::Stageable.
sub _needs_symlinks_followed_when_syncing {
    1;
}

1;
