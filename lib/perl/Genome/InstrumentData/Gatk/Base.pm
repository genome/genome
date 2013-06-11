package Genome::InstrumentData::Gatk::Base;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Gatk::Base {
    is => 'Genome::InstrumentData::AlignedBamResult',
    has_input => [
        bam_source => { # PROVIDES bam_file
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
        reference_build => { # PROVIDES fasta VIA full_consensus_path('fa')
            is => 'Genome::Model::Build::ImportedReferenceSequence',
        },
    ],
    has_constant => [
        _tmpdir => {  calculate => q| return File::Temp::tempdir(CLEANUP => 1); |, },
        # inputs
        input_bam_file => { 
            via => 'bam_source',
            to => 'bam_file', 
        },
        reference_fasta => { 
            calculate_from => [qw/ reference_build /],
            calculate => q| return $reference_build->full_consensus_path('fa'); |, 
        },
    ],
};

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $class = $self->class;
    $class =~ s/^Genome::InstrumentData::Gatk:://;
    $class =~ s/Result$//;
    return sprintf(
        "/alignment_data/gatk/%s-%s-%s-%s-%s", 
        Genome::Utility::Text::camel_case_to_string($class, '_'), 
        Sys::Hostname::hostname(),
        $ENV{USER}, $$, 1#$self->id,
    );
}

sub resolve_allocation_kilobytes_requested {
    my $self = shift;
    my $kb_requested = -s $self->input_bam_file;
    return int($kb_requested * 1.5);
}

sub resolve_allocation_disk_group_name {
    return 'info_alignments';
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $prepare_output_directory = eval{ $self->_prepare_output_directory; };
    if ( not $prepare_output_directory ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to prepare output directory!') if $@;
        return;
    }

    return $self;
}

1;

