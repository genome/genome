package Genome::Site::TGI::CaptureSet;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::CaptureSet {
    table_name => q|
        (select
            setup_name name,
            setup_id id,
            setup_status status,
            setup_description description
        from setup
        where setup_type = 'setup capture set'
        ) capture_set
    |,
    id_by => [
        id => { },
    ],
    has => {
        name => { },
        description => { },
        status => { },

        #Fields not native to CaptureSets to meet the interface of Genome::FeatureList::Command::Create:
        format => {
            is => 'Text',
            calculate => q{
                unless($self->_format) {
                    $self->_format($self->_resolve_format);
                }
                return $self->_format;
            },
            valid_values => Genome::FeatureList->__meta__->property('format')->valid_values,
        },
        source => {
            is => 'Text',
            calculate => q{
                unless($self->_source) {
                    $self->_source($self->_resolve_source);
                }
                return $self->_source;
            },
        },
        reference => {
            is => 'Number',
            calculate => q{
                unless($self->_reference) {
                    $self->_reference($self->_resolve_reference);
                }
                return $self->_reference;
            },
        },
        subject => {
            is => 'Number',
            calculate => q{
                return undef; #No way to determine at this time
            },
        },
        content_type => {
            is => 'Number',
            calculate => q{
                unless($self->_content_type) {
                    $self->_content_type($self->_resolve_content_type);
                }
                return $self->_content_type;
            },
        },
        file_path => {
            is => 'Text',
            calculate => q{
                unless($self->_file_path) {
                    $self->_file_path($self->_resolve_file_path);
                }
                return $self->_file_path;
            },
        },
    },
    has_optional => {
        file_storage_id => {
            calculate_from => ['_capture_set'],
            calculate => q{ $_capture_set->file_storage_id },
        },
    },
    has_transient => [
        _format => {
            is => 'Text',
            valid_values => Genome::FeatureList->__meta__->property('format')->valid_values,
        },
        _source => {
            is => 'Text',
        },
        _reference => {
            is => 'Genome::Model::Build::ReferenceSequence',
        },
        _file_path => {
            is => 'Text',
        },
        _content_type => {
            is => 'Text',
        },
    ],
    doc         => '',
    data_source => 'Genome::DataSource::Oltp',
};

sub _capture_set {
    return GSC::Setup::CaptureSet->get($_[0]->id);
}

sub barcodes {
    my $self = shift;
    my $cs = $self->_capture_set;
    my @barcodes = $cs->get_barcodes;
    return map {$_->barcode} @barcodes;
}

sub _resolve_format {
    my $self = shift;

    my $bed_file_content = $self->_capture_set->get_file_storage->content;

    unless($bed_file_content) {
        $self->error_message('Could not find BED file for capture set ' . $self->name . '.');
        die($self->error_message);
    }

    if($bed_file_content =~ /track name=/ims) {
        return 'multi-tracked';
    }

    my @lines = split("\n",$bed_file_content);
    for my $line (@lines) {
        chomp $line;
        my @fields = split("\t", $line);
        if($fields[1] == 0) {
            return 'true-BED';
        }
    }

    if($self->name =~ /agilent/i) {
        return '1-based';
    }

    if($self->name =~ /nimblegen/i) {
        return 'true-BED';
    }

    return 'unknown';
}

sub _resolve_source {
    my $self = shift;

    if($self->name =~ /agilent/i) {
        return 'agilent';
    } elsif ($self->name =~ /nimblegen/i) {
        return 'nimblegen';
    }

    return;
}

sub _resolve_file_path {
    my $self = shift;

    my $bed_file_content = $self->_capture_set->get_file_storage->content;

    unless($bed_file_content) {
        $self->error_message('Could not find BED file for capture set ' . $self->name . '.');
        die($self->error_message);
    }

    my $temp_bed_file = Genome::Sys->create_temp_file_path;
    Genome::Sys->write_file($temp_bed_file, $bed_file_content);

    return $temp_bed_file;
}

sub _resolve_reference {
    my $self = shift;

    #Try to determine reference based on original file name
    my $fs = $self->_capture_set->get_file_storage;
    if($fs) {
        my $name = $fs->file_name;
        if($name =~ /HG18/) {
            return Genome::Model::Build::ReferenceSequence->get_by_name('NCBI-human-build36');
        } elsif($name =~ /HG19/) {
            return Genome::Model::Build::ReferenceSequence->get_by_name('GRCh37-lite-build37');
        } else {
            return undef;
        }
    } else {
        return undef;
    }
}

sub _resolve_content_type {
    my $self = shift;

    return $self->_capture_set->content_type || undef;
}

1;
