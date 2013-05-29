package Genome::ProcessingProfile::AmpliconAssembly;

use strict;
use warnings;

use Genome;

use Bio::SeqIO;

class Genome::ProcessingProfile::AmpliconAssembly {
    is => 'Genome::ProcessingProfile::Staged',
    has_param => [
        sequencing_center => {
            doc => 'Place from whence the reads have come.',
            valid_values => [qw/ gsc broad baylor nisc /],
        },
        sequencing_platform => {
            doc => 'Platform (machine) from whence the reads where created.',
            valid_values => [qw/ sanger 454 solexa solid /],
        },
        assembler => {
            doc => 'Assembler type for assembling said reads.',
            valid_values => [qw/ maq newbler pcap phredphrap /],
        },
        assembly_size => {
            doc => 'Estimated assembly size, used for metrics and such',
        },
        region_of_interest => {
            doc => 'The name of the region being targeted',
        },
        purpose => { 
            doc => 'Purpose of these amplicon assemblies.',
            valid_values => [qw/ reference composition /],
        },
        primer_amp_forward => {
            doc => 'Primer used for amplification in the forward (5\') direction.  Enter both the name and sequence of the primer, separated by a colon as "NAME:SEQUENCE".  This will be used for naming the profile as well as orientation of the assemblies.', 
        },
        primer_amp_reverse => {
            doc => 'Primer used for amplification in the reverse (3\') direction.  Enter both the name and sequence of the primer, separated by a colon as "NAME:SEQUENCE".  This will be used for naming the profile as well as orientation of the assemblies.', 
        },
        primer_seq_forward => { 
            doc => 'Primer used for *internal* sequencing in the forward (5\') direction.  Enter both the name and sequence of the primer, separated by a colon as "NAME:SEQUENCE".  This will be used for naming the profile as well as orientation of the assemblies.', 
            is_optional => 1,
        },
        primer_seq_reverse => { 
            is_optional => 1,
            doc => 'Primer used for *internal* sequencing in the reverse (3\') direction.  Enter both the name and sequence of the primer, separated by a colon as "NAME:SEQUENCE".  This will be used for naming the profile as well as orientation of the assemblies.', 
        },
    ],
};

sub create {
    die __PACKAGE__ . ' is deprecated.';
}

sub delete {
    my $self = shift;

    # Remove primer files
    for my $method ( primer_fasta_methods() ) {
        unlink $self->$method if -e $self->$method;
    }

    return $self->SUPER::delete;
}

sub stages {
    return;
}

sub assemble_job_classes {
    return;
}

sub assemble_objects {
    return 1;
}

my %PRIMER_SENSES_AND_DIRECTIONS = (
    sense => 'forward',
    anti_sense => 'reverse',
);

sub primer_senses {
    return keys %PRIMER_SENSES_AND_DIRECTIONS;
}

sub primer_fasta_methods {
    return map { sprintf('%s_primer_fasta', $_) } primer_senses();
}

sub primer_fasta_directory {
    return '/gscmnt/839/info/medseq/processing_profile_data/amplicon_assembly';
}

sub sense_primer_fasta {
    return sprintf('%s/%s.sense.fasta', $_[0]->primer_fasta_directory, $_[0]->id);
}

sub anti_sense_primer_fasta {
    return sprintf('%s/%s.anti_sense.fasta', $_[0]->primer_fasta_directory, $_[0]->id);
}

1;
