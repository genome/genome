package Genome::Model::PhenotypeCorrelation::Command::DumpClinicalData;

use Genome;
use Sort::Naturally qw/nsort/;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::DumpClinicalData {
    is => "Command::V2",
    doc => "Dump clinical data for the specified samples (or population group) to a tsv file",
    has => [
        samples => {
            is => "Genome::Sample",
            doc => "List of samples (or a population group) by name or id",
            is_many => 1,
        },
        output_file => {
            is => "File",
            doc => "Path to the output file to write."
        },
        nomenclature => {
            is => "Genome::Nomenclature",
            doc => "nomenclature used to access clinical data",
        },
        missing_string => {
            is => "Text",
            doc => "What to output for missing fields (e.g., NA, null, ...)",
            default_value => "",
        }
    ],

    has_optional => [
        md5_file => {
            is => "File",
            doc => "If specified, the md5sum of the output file will be computed and stored in this file",
        },
    ],

    has_transient_optional => [
        _clin_fh => {
            is => 'IO::File',
        },
        _md5_fh => {
            is => 'IO::File',
        },
        _md5sum => {
            is => 'Text',
        },
    ],

};

sub help_synopsis {
    return <<"EOS"
    genome model phenotype-correlation dump-clinical-data \
        --samples <population group id> \
        --output-file /tmp/clindat.tsv \
        --nomenclature WUGC \
        --missing-string NA
EOS
}

sub help_detail {
    return <<"EOS"
Dump clinical attributes for the given samples to a file.
EOS
}

sub execute {
    my $self = shift;

    my @samples = $self->samples;
    my $n_samples = scalar(@samples);
    $self->status_message("Preparing clinical data file for $n_samples samples\n");
    my $clinical_data = Genome::Model::PhenotypeCorrelation::ClinicalData->from_database($self->nomenclature, @samples);

    # make sure we can open output files
    $self->_clin_fh(Genome::Sys->open_file_for_writing($self->output_file));


    $self->_md5sum($clinical_data->to_filehandle($self->_clin_fh, missing_string => $self->missing_string));
    $self->_clin_fh->close();
    $self->status_message("m5sum of clinical data: " . $self->_md5sum . "\n");

    if ($self->md5_file) {
        $self->_md5_fh(Genome::Sys->open_file_for_writing($self->md5_file));
        $self->_md5_fh->write($self->_md5sum . "\n");
        $self->_md5_fh->close();
    }

    return 1;
}

1;
