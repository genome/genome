package Genome::Model::Build::ReferenceSequence::AlignerIndex;

use Genome;
use warnings;
use strict;


class Genome::Model::Build::ReferenceSequence::AlignerIndex {
    is => ['Genome::SoftwareResult::Stageable'],

    has => [

        reference_build         => {
                                    is => 'Genome::Model::Build::ImportedReferenceSequence',
                                    id_by => 'reference_build_id',
                                },
        reference_name          => { via => 'reference_build', to => 'name', is_mutable => 0, is_optional => 1 },

        aligner                 => {
                                    calculate_from => [qw/aligner_name aligner_version aligner_params/],
                                    calculate => q|no warnings; "$aligner_name $aligner_version $aligner_params"|
                                },
    ],
    has_input => [
        reference_build_id      => {
                                    is => 'Number',
                                    doc => 'the reference to use by id',
                                },
    ],
    has_param => [
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
                                    doc => 'any additional params for the aligner in a single string',
                                },
    ],

    has_transient => [
        aligner_class_name      => {
                                    is => 'Text',
                                    is_optional => 1,
        }
    ]
};

sub _working_dir_prefix {
    "aligner-index";
}

sub __display_name__ {
    my $self = shift;
    my @class_name = split("::", $self->class);
    my $class_name = $class_name[-1];
    no warnings;
    return sprintf("%s for build %s with %s, version %s, params='%s'",
        $class_name,
        $self->reference_name,
        $self->aligner_name,
        $self->aligner_version,
        $self->aligner_params || "");
}

sub required_rusage {
    # override in subclasses
    # e.x.: "-R 'select[model!=Opteron250 && type==LINUX64] span[hosts=1] rusage[tmp=50000:mem=12000]' -M 1610612736";
    ''
}

sub aligner_requires_param_masking {
    my $class = shift;
    my $aligner_name = shift;

    my $aligner_class = 'Genome::InstrumentData::AlignmentResult::'  . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($aligner_name);

    # if aligner params are not required for index, and we can   generically create an index for that version, then filter it out.
    if ($aligner_class->aligner_params_required_for_index) {
        $class->status_message("This aligner requires a parameter-specific index.  Can't mask params out.");
        return 0;
    }

    return 1;
}

sub _supports_multiple_reference {
    my $self = shift;
    my $aligner_name = $self->aligner_name;
    my $aligner_class = 'Genome::Model::Tools::'  . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($aligner_name);
    return unless $aligner_class->can('supports_multiple_reference');
    return $aligner_class->supports_multiple_reference($self->aligner_version);
}

sub get_or_create {
    my $class = shift;

    my @objects = $class->SUPER::get_or_create(@_);

    for my $obj (@objects) {
        next unless ref($obj); # sometimes UR gives us back the package name when deleting?
        unless ($obj->generate_dependencies_as_needed()) {
            $obj->error_message("Failed to get AlignmentIndex objects for dependencies of " . $obj->__display_name__);
            return;
        }
    }

    if (@objects > 1) {
        return @objects if wantarray;
        my @ids = map { $_->id } @objects;
        die "Multiple matches for $class but get or create was called in scalar context! Found ids: @ids";
    }
    else {
        return $objects[0];
    }
}

sub create {
    my $class = shift;
    my %p = @_;

    my $aligner_class = 'Genome::InstrumentData::AlignmentResult::'  . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($p{aligner_name});
    $class->status_message(sprintf("Resolved aligner class %s, making sure it's real and can be loaded.", $aligner_class));
    unless ($aligner_class->class) {
        $class->error_message(sprintf("Failed to load aligner class (%s).", $aligner_class));
        return;
    }

    my $self = $class->SUPER::create(%p);
    return unless $self;
    $self->aligner_class_name($aligner_class);

    $self->status_message("Prepare staging directories...");
    unless ($self->_prepare_staging_directory) {
        $self->error_message("Failed to prepare working directory");
        return;
    }

    unless ($self->_prepare_reference_index) {
        $self->error_message("Failed to prepare reference index!");
        return;
    }

    unless ($self->generate_dependencies_as_needed()) {
        $self->error_message("Failed to create AlignmentIndex objects for dependencies");
        return;
    }

    return $self;
}

# TODO:
# get() calls this method, and has a side-effect of creating dependent aligner indexes
# 1. if you have side effects (avoid in general where possible), don't put them in a method called check_*
# 2. don't override get(), make another method with the combined effect of getting data and doing work
# -ssmith
sub generate_dependencies_as_needed {
    my $self = shift;

    # if the reference is a compound reference
    if ($self->reference_build->append_to) {
        my %params = (
            aligner_name => $self->aligner_name,
            aligner_params => $self->aligner_params,
            aligner_version => $self->aligner_version,
        );

        for my $b ($self->reference_build->append_to) { # (append_to is_many)
            $params{reference_build} = $b;
            $self->status_message("Creating AlignmentIndex for build dependency " . $b->name);
            my $result = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_or_create(%params);
            unless($result) {
                die $self->error_message("Failed to create AlignmentIndex for dependency " . $b->name);

            }
        }
    }

    return 1;
}

sub _prepare_reference_index {
    my $self = shift;

    my $reference_fasta_file;
    if ($self->_supports_multiple_reference) {
        $reference_fasta_file = $self->reference_build->primary_consensus_path('fa');
    } else {
        $reference_fasta_file = $self->reference_build->full_consensus_path('fa');
    }

    unless (-s $reference_fasta_file) {
        $self->error_message(sprintf("Reference fasta file %s does not exist", $reference_fasta_file));
        return;
    }

    $self->status_message(sprintf("Confirmed non-zero reference fasta file is %s", $reference_fasta_file));
    unless (symlink($reference_fasta_file, sprintf("%s/all_sequences.fa", $self->temp_staging_directory))) {
        $self->error_message("Couldn't symlink reference fasta into the staging directory");
    }

    my $reference_remap_file = sprintf("%s.remap", $reference_fasta_file);

    if (-e $reference_remap_file) {
        $self->status_message("Detected $reference_remap_file.remap. Symlinking that as well.");
        unless (symlink($reference_remap_file, sprintf("%s/all_sequences.fa.remap", $self->temp_staging_directory))) {
            $self->error_message("Couldn't symlink reference remap into the staging directory");
        }
    }

    unless ($self->aligner_class_name->prepare_reference_sequence_index($self)) {
        $self->error_message("Failed to prepare reference sequence index.");
        return;
    }

    my $output_dir = $self->output_dir || $self->_prepare_output_directory;
    $self->status_message("Alignment output path is $output_dir");

    unless ($self->_promote_data)  {
        $self->error_message("Failed to de-stage data into output path " . $self->output_dir);
        return;
    }

    $self->_reallocate_disk_allocation;

    $self->status_message("Prepared alignment reference index!");

    return $self;
}

sub _modify_params_for_lookup_hash {
    my ($class, $params_ref) = @_;

    if (exists $params_ref->{aligner_name} &&
            $class->aligner_requires_param_masking($params_ref->{aligner_name})) {
        $params_ref->{aligner_params} = undef;
    }
}

sub _gather_params_for_get_or_create {
    my $class = shift;
    my $p = $class->SUPER::_gather_params_for_get_or_create(@_);

    unless ($p->{params}{test_name}) {
        $p->{params}{test_name} = ($ENV{GENOME_ALIGNER_INDEX_TEST_NAME} || undef);
    }
    if (exists $p->{params}{aligner_name} && $class->aligner_requires_param_masking($p->{params}{aligner_name})) {
        $p->{params}{aligner_params} = undef;
    }

    #my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %{ $p->{inputs} });
    #my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %{ $p->{params} });
    #$p->{software_result_params}{params_id} = $params_bx->id;
    #$p->{software_result_params}{inputs_id} = $inputs_bx->id;

    return $p;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $aligner_name_tag = $self->aligner_name;
    $aligner_name_tag =~ s/[^\w]/_/g;

    my $staged_basename = File::Basename::basename($self->temp_staging_directory);
    my @path_components = ('model_data','ref_build_aligner_index_data',$self->reference_build->model->id,'build'.$self->reference_build->id, $staged_basename);

    if ($self->test_name) {
        push @path_components, "test_".$self->test_name;
    }

    push @path_components, $aligner_name_tag;

    my $aligner_version_tag = $self->aligner_version;
    $aligner_version_tag =~ s/[^\w]/_/g;
    push @path_components, $aligner_version_tag;

    if ($self->aligner_params) {
        my $aligner_params_tag = $self->aligner_params;
        $aligner_params_tag =~ s/[^\w]/_/g;
        push @path_components, $aligner_params_tag;
    }

    my $directory = join('/', @path_components);

    $self->status_message(sprintf("Resolved allocation subdirectory to %s", $directory));
    return $directory;
}

sub resolve_allocation_disk_group_name {
    if ($_[0]->reference_build->model->is_rederivable) {
        return 'info_genome_models';
    } else {
        return $ENV{GENOME_DISK_GROUP_REFERENCES};
    }
}

sub full_consensus_path {
    my $self = shift;
    my $extension = shift;

    return $self->output_dir . "/all_sequences." . $extension;
}




1;
