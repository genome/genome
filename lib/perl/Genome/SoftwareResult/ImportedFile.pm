package Genome::SoftwareResult::ImportedFile;

use strict;
use warnings;
use Genome;

class Genome::SoftwareResult::ImportedFile {
    is => 'Genome::SoftwareResult::StageableSimple::SingleFile',
    has_input => [
        file_content_hash => {
            is => 'Text',
            doc => 'MD5 hash of the file that is contained in this result',
        },
    ],
    has_optional_metric => [
        source_file_path => {
            is => 'Path',
            doc => 'Path to the original file used to create this result',
        }
    ],
};

sub _modify_params_for_lookup_hash {
    my ($class, $params_ref) = @_;

    unless (defined $params_ref->{file_content_hash}) {
        if (my $file = delete $params_ref->{source_file_path}) {
            $file = Cwd::abs_path($file);
            $class->fatal_message('Source file path (%s) given to create software result does not exist!', $file) if not -s $file;

            my $md5 = Genome::Sys->md5sum($file);
            $class->fatal_message('No md5 for source file! %s', $file) if not $md5;

            $params_ref->{'file_content_hash'} = $md5;
        }
    }
}

sub _gather_params_for_get_or_create {
    my $class = shift;

    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);
    if($bx->specifies_value_for('source_file_path')) {
        my $file = $bx->value_for('source_file_path');
        $bx = $bx->add_filter('file_content_hash', Genome::Sys->md5sum($file));
    }

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

    return {
        inputs => \%is_input,
        params => \%is_param,
    };
}

sub _run{
    my $self = shift;
    $self->file_content_hash(Genome::Sys->md5sum($self->source_file_path));
    $self->lookup_hash($self->calculate_lookup_hash);
    Genome::Sys->create_symlink(Cwd::abs_path($self->source_file_path), $self->_temp_staging_file_path);
}

sub _file_name {
    my $self = shift;
    return File::Basename::basename($self->source_file_path);
}

sub _needs_symlinks_followed_when_syncing { return 1; }

1;
