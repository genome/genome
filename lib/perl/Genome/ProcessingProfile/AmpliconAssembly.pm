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
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;


    # Check primers, create sense fastas and create name
    for my $sense ( primer_senses() ) {
        my $file_method = sprintf('%s_primer_fasta', $sense);
        my $file = $self->$file_method;
        unlink $file if -e $file;
        my $bioseq_io = Bio::SeqIO->new(
            '-file' => ">$file",
            '-format' => 'Fasta',
        );
        for my $primer_purpose ( primer_purposes() ) {
            my $method = sprintf('primer_%s_%s', $primer_purpose, _primer_direction_for_sense($sense));
            my $primer = $self->$method;
            next unless defined $primer;
            my ($primer_name, $primer_seq) = split(/:/, $primer);
            unless ( $primer =~ /^([\w\d\-\.\_]+):([ATCGMRWSYKVHDBN]+)$/i ) {
                $self->error_message("Invlaid format for primer ($primer).  Please use 'NAME:SEQUENCE'");
                $self->delete;
                return;
            }
            $bioseq_io->write_seq( 
                Bio::PrimarySeq->new( 
                    '-id' => $1,
                    '-seq' => $2,
                    '-alphabet' => 'dna',
                )
            );
        }
    }

    unless ( grep { defined $self->$_ and -s $self->$_ } primer_fasta_methods() ) {
        $self->error_message("No primers indicated for processing profile");
        $self->delete;
        return;
    }

    return $self;
}

sub delete {
    my $self = shift;

    # Remove primer files
    for my $method ( primer_fasta_methods() ) {
        unlink $self->$method if -e $self->$method;
    }

    return $self->SUPER::delete;
}

#< BUILDING >#
sub stages {
    return (qw/
        assemble
        /);
}

sub assemble_job_classes {
    return (qw/
        Genome::Model::Event::Build::AmpliconAssembly::VerifyInstrumentData
        Genome::Model::Event::Build::AmpliconAssembly::PrepareInstrumentData
        Genome::Model::Event::Build::AmpliconAssembly::TrimAndScreen
        Genome::Model::Event::Build::AmpliconAssembly::Assemble
        Genome::Model::Event::Build::AmpliconAssembly::Classify
        Genome::Model::Event::Build::AmpliconAssembly::Orient
        Genome::Model::Event::Build::AmpliconAssembly::Collate
        Genome::Model::Event::Build::AmpliconAssembly::CleanUp
        Genome::Model::Event::Build::AmpliconAssembly::Reports
        /);
}

sub assemble_objects {
    return 1;
}

#< Primers >#
my %PRIMER_SENSES_AND_DIRECTIONS = (
    sense => 'forward',
    anti_sense => 'reverse',
);
sub primer_senses {
    return keys %PRIMER_SENSES_AND_DIRECTIONS;
}

sub primer_directions {
    return values %PRIMER_SENSES_AND_DIRECTIONS;
}

sub _primer_direction_for_sense {
    die "Need primer sense to get direction\n" if not defined $_[0] or $_[0] eq __PACKAGE__;
    return $PRIMER_SENSES_AND_DIRECTIONS{$_[0]};
}

sub primer_purposes {
    return (qw/ amp seq /);
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

#$HeadURL$
#$Id$
