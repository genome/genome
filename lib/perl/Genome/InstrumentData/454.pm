package Genome::InstrumentData::454;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::454 {
    is => 'Genome::InstrumentData',
    has_constant => [
        sequencing_platform => { value => '454' },
    ],
    has_optional => [
        sff_file => { is_attribute => 1, },
        beads_loaded => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'beads_loaded' ],
            is_mutable => 1,
        },
        copies_per_bead => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'copies_per_bead' ],
            is_mutable => 1,
        },
        fc_id => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'fc_id' ],
            is_mutable => 1,
        },
        incoming_dna_name => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'incoming_dna_name' ],
            is_mutable => 1,
        },
        key_pass_wells => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'key_pass_wells' ],
            is_mutable => 1,
        },
        predicted_recovery_beads => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'predicted_recovery_beads' ],
            is_mutable => 1,
        },
        region_id => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'region_id' ],
            is_mutable => 1,
        },
        region_number => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'region_number' ],
            is_mutable => 1,
        },
        research_project => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'research_project' ],
            is_mutable => 1,
        },
        sample_set => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'sample_set' ],
            is_mutable => 1,
        },
        ss_id => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'ss_id' ],
            is_mutable => 1,
        },
        supernatant_beads => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'supernatant_beads' ],
            is_mutable => 1,
        },
        total_key_pass => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'total_key_pass' ],
            is_mutable => 1,
        },
        total_raw_wells => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'total_raw_wells' ],
            is_mutable => 1,
        },
        index_sequence => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'index_sequence' ],
            is_mutable => 1,
        },
        read_count => { # the universal way to get number of reads
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'total_reads' ], # doesn't match ok
            is_mutable => 1,
        },
        total_reads => { # old method
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'total_reads' ],
            is_mutable => 1,
        },
        total_bases_read => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'total_bases_read' ],
            is_mutable => 1,
        },
        is_paired_end => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'paired_end' ],
            is_mutable => 1,
        },
        ],
        has_optional_transient => [
        _fasta_file => { is => 'FilePath', is_mutable => 1, },
        _qual_file => { is => 'FilePath', is_mutable => 1, },
        ],
    };

    BEGIN: {
        Genome::InstrumentData::Solexa->class;
        no warnings 'once';
        *dump_trimmed_fastq_files = \&Genome::InstrumentData::Solexa::dump_trimmed_fastq_files;
    }

    sub full_path {
    Carp::confess("Full path is not valid for 454 instrument data");
}

sub bam_path {
    my $self = shift;
    $self->warning_message("Asked this 454 instrument data for bam path, but this is not implemented yet.");
    return undef;
}

sub _default_full_path {
    my $self = shift;
    return sprintf('%s/%s/%s', $self->_data_base_path, $self->run_name, $self->region_id);
}

sub calculate_alignment_estimated_kb_usage {
    my $self = shift;
    return 500000;
}

sub is_external {
    return;
}

#< Fastq, Fasta, Qual ... >#
sub dump_sanger_fastq_files {
    my $self = shift;
    
    my %params = @_;
    
    unless (-s $self->sff_file) {
        $self->error_message(sprintf("SFF file (%s) doesn't exist for 454 instrument data %s", $self->sff_file, $self->id));
        die $self->error_message;
    }
    
    my $dump_directory = delete $params{'directory'} || Genome::Sys->base_temp_directory();
    
    my $output_file = sprintf("%s/%s-output.fastq", $dump_directory, $self->id);
    
    my %sff2fastq_params = (
        sff_file => $self->sff_file,
        fastq_file => $output_file,
    );
    $sff2fastq_params{force_untrimmed} = $params{force_untrimmed} if exists $params{force_untrimmed};
    my $cmd = Genome::Model::Tools::454::Sff2Fastq->create(%sff2fastq_params);
    unless ($cmd->execute) {
        $self->error_message("Sff2Fastq failed while dumping fastq file for instrument data " . $self->id);
        die $self->error_message;
    }
    
    unless (-s $output_file) {
        $self->error_message("Sff2Fastq claims it worked, but the output file was gone or empty length while dumping fastq file for instrument data "
                             . $self->id . " expected output file was $output_file");
        die $self->error_message;
    }
    
    return ($output_file);
}

sub dump_fasta_file {
    my ($self, %params) = @_;
    $params{type} = 'fasta';
    return $self->_run_sffinfo(%params);
}

sub dump_qual_file {
    my ($self, %params) = @_;
    $params{type} = 'qual';
    return $self->_run_sffinfo(%params);
}

sub _run_sffinfo {
    my ($self, %params) = @_;

    # Type 
    my $type = delete $params{type};
    my %types_params = (
        fasta => '-s',
        qual => '-q',
    );
    unless ( defined $type and grep { $type eq $_ } keys %types_params ) { # should not happen
        Carp::confess("No or invalid type (".($type || '').") to run sff info.");
    }

    # Verify 64 bit
    unless ( Genome::Config->arch_os =~ /x86_64/ ) {
        Carp::confess(
            $self->error_message('Dumping $type file must be run on 64 bit machine.')
        );
    }
    
    # SFF
    my $sff_file = $self->sff_file;
    unless ( -s $sff_file ) {
        Carp::confess(
            $self->error_message(
                "SFF file ($sff_file) doesn't exist for 454 instrument data (".$self->id.")"
            )
        );
    }

    # File
    my $directory = delete $params{'directory'} 
        || Genome::Sys->base_temp_directory();
    my $file = sprintf("%s/%s.%s", $directory, $self->id, $type);
    unlink $file if -e $file;
    
    # SFF Info
    my $sffinfo = Genome::Model::Tools::454::Sffinfo->create(
        sff_file => $sff_file,
        output_file => $file,
        params => $types_params{$type},
    );
    unless ( $sffinfo ) {
        Carp::confess(
            $self->error_message("Can't create SFF Info command.")
        );
    }
    unless ( $sffinfo->execute ) {
        Carp::confess(
            $self->error_message("SFF Info command failed to dump $type file for instrument data (".$self->id.")")
        );
    }

    # Verify
    unless ( -s $file ) {
        Carp::confess(
            $self->error_message("SFF info executed, but a fasta was not produced for instrument data (".$self->id.")")
        );
    }

    return $file;
}

#< Run Info >#
sub run_identifier {
    my $self = shift;
    return sprintf('%s.%s%s', $self->run_name, $self->region_number, ( $self->index_sequence ? '-'.$self->index_sequence : '' ));
}

sub run_start_date_formatted {
    my $self = shift;

    my ($y, $m, $d) = $self->run_name =~ m/R_(\d{4})_(\d{2})_(\d{2})/;

    my $dt_format = UR::Time->config('datetime');
    UR::Time->config(datetime=>'%Y-%m-%d');
    my $dt = UR::Time->numbers_to_datetime(0, 0, 0, $d, $m, $y);
    UR::Time->config(datetime=>$dt_format);

    return $dt; 
}
#<>#

1;

