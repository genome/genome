package Genome::Sample::Command::Import::Manager;

use strict;
use warnings;

use Genome;

use Data::Dumper;
use IO::File;

class Genome::Sample::Command::Import::Manager {
    is => 'Command::V2',
    doc => 'Manage importing a group of samples including importing instrument data and creating and building models.',
    has_optional => [
        working_directory => {
            doc => 'Directory to read and write.',
        },
    ],
    has_calculated => [
        sample_csv_file => {
            calculate_from => 'working_directory',
            calculate => sub{ my $working_directory = shift; return $working_directory.'/samples.csv'; },
            doc => 'CSV file of samples and attributes. A column called "name" is required. The name should be dash (-) separated values of the nomenclature, indivdual id and sample id.',
        },
        job_dispatch_config => {
            calculate_from => 'working_directory',
            calculate => sub{ my $working_directory = shift; return $working_directory.'/job_dispatch.config'; },
        },
    ],
};

sub execute {
    my $self = shift;

    my $samples = $self->_load_samples;
    return if not $samples;

    $self->status($samples);

    return 1;
}

sub _load_samples {
    my $self = shift;

    my $samples = $self->_load_samples_from_csv_file;
    return if not $samples;

    my $load_job_status = $self->_load_job_status($samples);
    return if not $load_job_status;

    for my $sample ( 
        Genome::Sample->get(
            'name in' => [ keys %$samples ],
            '-hint' => [qw/ models instrument_data /],
        )
    ) {
        my $name = $sample->name;
        $samples->{$name}->{sample} = $sample;
        $samples->{$name}->{id} = $sample->id;
        $samples->{$name}->{inst_data} = ($sample->instrument_data)[-1];
        $samples->{$name}->{bam_path} = eval{ $sample->{inst_data}->bam_path };
        $samples->{$name}->{model} = ($sample->models)[-1];
        $samples->{$name}->{build} = eval{ $sample->{model}->latest_build };
    }

    return [ sort { $a->{name} cmp $b->{name} } values %$samples ];
}

sub _load_samples_from_csv_file {
    my $self = shift;
    
    my $sample_csv_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->sample_csv_file,
        separator => ',',
    );
    if ( not $sample_csv_reader ) {
        $self->error_message('Failed to open sample csv! '.$self->sample_csv_file);
        return;
    }

    if ( not grep { $_ eq 'name' } @{$sample_csv_reader->headers} ) {
        $self->error_message('No name column in sample csv! '.$self->sample_csv_file);
        return;
    }


    my %samples;
    while ( my $sample = $sample_csv_reader->next ) {
        my $name = delete $sample->{name};
        $samples{$name}->{name} = $name;
        $samples{$name}->{from_csv} = $sample;
    }

    return \%samples;
}

sub _load_job_status {
    my ($self, $samples) = shift;

    Carp::confess('Need samples to load job status!') if not S$samples;

    my $fh = IO::File->new('bjobs -g/ebelter/whisp_imports -w 2> /dev/null |');
    $fh->getline;
    while ( my $line = $fh->getline ) {
        my @tokens = split(/\s+/, $line);
        next if not $samples->{$tokens[6]};
        $samples->{$tokens[6]}->{job_status} = lc $tokens[2];
    }

    return 1;
}

sub status {
    my ($self, $samples) = @_;
    $self->status_message('Status');
    return $self->_status($samples);
}

sub _status {
    my %totals;
    my $status;
    for my $sample ( sort { $a->{name} cmp $b->{name} } @_ ) {
        $totals{total}++;
        $sample->{status} = _get_status($sample);
        $totals{ $sample->{status} }++;
        $totals{build}++ if $sample->{status} =~ /^build/;
        $status .= sprintf("%-20s %10s\n", $sample->{name}, $sample->{status});
    }
    print "$status\nSummary:\n".join("\n", map { sprintf('%-16s %s', $_, $totals{$_}) } sort { $a cmp $b } keys %totals)."\n";
    return 1;
}

sub cleanup_failed {
    my @samples = _load_samples();
    for my $sample ( @samples ) {
        $sample->{status} = _get_status($sample);
        next if $sample->{status} ne 'import_failed' or ( defined $sample->{bam_path} and -s $sample->{bam_path} );
        print $sample->{name}.' '.$sample->{status}." REMOVE!\n";
        $sample->{inst_data}->delete;
    }
    UR::Context->commit;
    return 1;
}

sub commands {
    print STDERR "Commands...\n";
    my @samples = _load_samples();
    for my $sample ( @samples ) {
        $sample->{status} = _get_status($sample);
        next if not _needs_import($sample);
        my $cmd = 'grep '.$sample->{name}.' sra_import_commands';
        print `$cmd`;
    }
    return 1;
}

sub models {
    print STDERR "Models...\n";
    my @samples = _load_samples();
    my %params = (
        processing_profile_id => 2673537,
        reference_sequence_build => Genome::Model::Build->get(106942997),
        dbsnp_build => Genome::Model::Build->get(106375969),
    );
    for my $sample ( @samples ) {
        if ( $sample->{model} ) {
            if ( $sample->{inst_data} and not $sample->{model}->instrument_data ) {
                $sample->{model}->add_instrument_data($sample->{inst_data});
            }
            if ( $sample->{build} and lc($sample->{build}->status) eq 'unstartable' ) {
                $sample->{model}->build_requested(1);
            }
            elsif ( not $sample->{build} and not $sample->{model}->build_requested ) {
                $sample->{model}->build_requested(1);
            }
            elsif ( $sample->{build} and $sample->{model}->build_requested ) {
                $sample->{model}->build_requested(0);
            }
            next;
        }
        $sample->{model} = Genome::Model->create(
            subject_id => $sample->{id},
            %params,
        );
        die 'Failed to create model! '.$sample->{name} if not $sample->{model};
        $sample->{model}->build_requested(1);
    }
    UR::Context->commit;
    return _status(@samples);
}

sub _needs_import {
    my $sample = shift;
    return 1 if $sample->{status} eq 'import_failed' or $sample->{status} eq 'import_needed';
}

sub _get_status {
    my $sample = shift;
    return 'import_'.$sample->{job_status} if $sample->{job_status};
    return 'import_needed' if not $sample->{inst_data};
    return 'import_failed' if not defined $sample->{bam_path} or not -s $sample->{bam_path};
    return 'model_needed' if $sample->{inst_data} and not $sample->{model};
    return 'build_requested' if $sample->{model}->build_requested;
    return 'build_needed' if not $sample->{build};
    return 'build_'.lc($sample->{build}->status);
}

# Imported bam and SRA, successful build for bam: sample.name=dbGaP-295277-661565
