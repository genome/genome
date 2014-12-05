package Genome::VariantReporting::Command::Wrappers::CreateIgvSession;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Framework::FileLookup qw(calculate_lookup);
use Digest::MD5 qw(md5_hex);

my $_JSON_CODEC = new JSON->allow_nonref;

class Genome::VariantReporting::Command::Wrappers::CreateIgvSession {
    is => 'Genome::Command::DelegatesToResult',
    has_input => [
        bam_hash_json => {
            is => 'String',
            doc => 'Hash where keys are the labels for the bam track and the values are the bam file paths',
        },
        genome_name => {
            is => 'String',
            doc => 'Name for the file and tracks',
        },
        merged_bed_reports => {
            is => 'Genome::VariantReporting::Framework::Component::Report::MergedReport',
            is_many => 1,
            doc => 'Reports to get the bed files from',
        },
        reference_name => {
            is => 'String',
            doc => 'name of the reference sequence build',
        },
        user => {
            is => 'Genome::Process',
            is_optional => 1,
            id_by => 'process_id',
        },
        process_id => {
            is => 'Text',
        },
    ],
};

sub result_class {
    return "Genome::VariantReporting::Command::Wrappers::IgvSession";
}

sub input_hash {
    my $self = shift;
    return (
        bam_hash_json => $self->bam_hash_json,
        bam_hash_json_lookup => $self->calculate_bam_hash_lookup,
        genome_name => $self->genome_name,
        merged_bed_reports => [$self->merged_bed_reports],
        reference_name => $self->reference_name,
    );
}

sub calculate_bam_hash_lookup {
    my $self = shift;
    my %bam_hash = %{$_JSON_CODEC->decode($self->bam_hash_json)};

    my %encoded_path_hash;

    while (my ($label, $bam_path) = each %bam_hash) {
        $encoded_path_hash{$label} = calculate_lookup($bam_path);
    }
    return md5_hex($_JSON_CODEC->canonical->encode(%encoded_path_hash));
}

1;

