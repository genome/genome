package Genome::Model::Tools::FastTier::MakeTierBitmasks;

use strict;
use warnings;
use Bit::Vector;
use Genome;
use UR;
use IO::File;
use Data::Dumper;

class Genome::Model::Tools::FastTier::MakeTierBitmasks {
    is => 'Command',
    has => [
        output_directory => {
            type => 'Text',
            is_input => 1,
            doc => 'Location for tier bitmasks to be dropped',
        },
        reference_sequence_build => {
            type => 'Genome::Model::Build::ReferenceSequence',
            is_input => 1,
            doc => 'Reference sequence to use for tier mask creation, default is NCBI human build36',
        },
        transcript_version => {
            type => 'Text',
            is_input => 1,
            doc => 'which version of the transcript to use',
        },
        ucsc_directory => {
            type => 'Text',
            is_input => 1,
            doc => 'The location of phastcons17,28, regulatory regions, etc',
        },
        species => {
            type => 'Text',
            is_input => 1,
        },
        build => {
            is_input => 1,
        },
        annotation_import_version => {
            is_input => 1,
        },
    ],
};


sub help_brief {
    "Used to generate bitmask files of teirs 1 through 4"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools fast-tier make-tier-bitmask ...    
EOS
}

sub execute {
    my $self = shift;

    my $build = $self->build;
    my $users = Genome::SoftwareResult::User->user_hash_for_build($build);
    $users->{'fast-tier tier-bitmasks'} = $build;

    my $result = Genome::Model::Tools::FastTier::TierBitmasks->get_or_create(
        reference_sequence_build => $self->reference_sequence_build,
        annotation_structures => Genome::Db::Ensembl::AnnotationStructures->get_or_create(version => $self->transcript_version,
               software_version => $self->annotation_import_version,
               reference_build_id => $self->reference_sequence_build->id,
               species => $self->species,
               data_set => 'Core',
               test_name => (Genome::Config::get('software_result_test_name') || undef),
        ),
        ucsc_directory => $self->ucsc_directory,
        test_name => (Genome::Config::get('software_result_test_name') || undef),
        species => $self->species,
        user_data => $users,
    );

    unless ($result) {
        $self->error_message("Failed to generate TierBitmasks result");
        return;
    }

    $self->debug_message("Using TierBitmasks result: ".$result->id);

    foreach my $path ($result->result_paths) {
        Genome::Sys->create_symlink(
            $result->output_dir."/".$path,
            $self->output_directory."/".$path
        );
    }

    return 1;
}

