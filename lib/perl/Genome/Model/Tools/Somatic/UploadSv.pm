package Genome::Model::Tools::Somatic::UploadSv;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Somatic::UploadSv {
    is => 'Command',
    has => [
    breakdancer_file => {
        is  => 'String',
        doc => 'The file of somatic pipeline results to be uploaded. This will usually be a high confidence tier 1 or 2 snp file, or a tier 1 indel file from the somatic pipeline.',
    },
    build_id => {
        is => 'Number',
        doc => 'The build id that should be linked to the variant. This is manual for now and required.',
    },
    no_header => {
        type => 'Boolean',
        default => 0,
        doc => 'If there is no header line in the file pass in a 1 and the first line will not be skipped.',
    },
    ],
};

sub help_brief {
    "Adds structural variants from breakdancer to the database tables for tracking SV" 
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    gmt somatic upload-sv --breakdancer-file sv.out --build-id 12345 
EOS
}

sub help_detail {                           
    return <<EOS
Adds structural variants from breakdancer to the database tables for tracking SV
EOS
}

sub execute {
    my $self = shift;

    my $variant_fh = IO::File->new($self->breakdancer_file);
    unless ($variant_fh) {
        $self->error_message("Could not open variant file: " . $self->breakdancer_file . " for reading. $!");
        die;
    }

    # Toss the header if it is present
    unless ($self->no_header) {
        my $header = $variant_fh->getline;
    }

    # For each SV in the file, make a new SV object if one does not already exist. If one already exists it should be exactly the same so we will just make a new / updated build link to it later
    while (my $line = $variant_fh->getline) {
        chomp $line;
        # Currently "run_param" seems to be empty... will upload as null to the DB but watch out for future column additions 
        my ($chromosome_1, $position_1, $orientation_1, $chromosome_2, $position_2, $orientation_2, $event_type, $sv_size, $score, $num_reads, $num_reads_lib, $allele_frequency, $version, $run_param) = split ("\t", $line);

        # Create or grab the existing SV to add a build link later
        my $new_sv = Genome::Model::SV->get_or_create(
            chromosome_1    => $chromosome_1,
            pos_predicted_1 => $position_1,
            orientation_1   => $orientation_1,
            pos_predicted_2 => $position_2,
            chromosome_2    => $chromosome_2,
            orientation_2   => $orientation_2,
            event_type      => $event_type,
            sv_size         => $sv_size,
        );
        unless($new_sv) {
            $self->error_message("UR was unable to instantiate a variant object for this line:");
            $self->error_message($line);
            die;
        }

        my $build = Genome::Model::Build->get($self->build_id);
        unless ($build) {
            $self->error_message("No build found for build_id " . $self->build_id . "... please supply a valid build id");
            die;
        }

        # Get the old SV to build link and delete it, since we are creating a new one with updated information
        # This will only matter if we are rerunning a build which we probably shouldnt do anyways.
        my $existing_build_sv = Genome::Model::BuildSV->get(
            sv                 => $new_sv,
            build              => $build,
        );
        if ($existing_build_sv) {
            $existing_build_sv->delete;
        }

        my $new_build_sv = Genome::Model::BuildSV->create(
            sv                 => $new_sv,
            build              => $build,
            allele_frequency   => $allele_frequency,
            breakdancer_score  => $score,
            num_reads          => $num_reads,
            num_reads_lib      => $num_reads_lib,
            run_param          => $run_param,
            version            => $version,
            somatic_status     => 'somatic', # FIXME hard code this for now... obviously will not always be the case
        );
        
        # TODO add sv validation information to mirror the snp/indel upload tool once the table and class has been made...
        # Ken Chen is still working on this format for now

        unless($new_build_sv) {
            $self->error_message("Unable to create Build-Variant Link");
            $self->error_message("Problem line: $line");
            die;
        }
    }

    return 1;
}

1;
