package Genome::Sample::Command::Import::Manager;

use strict;
use warnings;

use Genome;

use Data::Dumper;
use IO::File;

class Genome::Sample::Command::Import::Manager {
    is => 'Command::V2',
    doc => 'Manage importing a group of samples including importing instrument data and creating and building models.',
    has => [
    ],
};

sub execute {
    my $self = shift;

    my $samples = $self->_load_samples;
    return if not $samples;

    return 1;
}

sub _load_samples {
    my $self = shift;

    my $fh = IO::File->new('bjobs -g/ebelter/whisp_imports -w 2> /dev/null |');
    $fh->getline;
    my %jobs_status;
    while ( my $line = $fh->getline ) {
        my @tokens = split(/\s+/, $line);
        $jobs_status{$tokens[6]} = lc $tokens[2];
    }

    my @samples;
    for my $sample ( 
        Genome::Sample->get(
            'name in' => [ _sample_names() ],
            '-hint' => [qw/ models instrument_data /],
        )
    ) {
        my $inst_data = ($sample->instrument_data)[-1];
        my $bam_path = eval{ $inst_data->bam_path };
        my $model = ($sample->models)[-1];
        my $build = eval{ $model->latest_build };
        push @samples, {
            name => $sample->name,
            id => $sample->id,
            sample => $sample,
            job_status => $jobs_status{$sample->name},
            inst_data => $inst_data,
            bam_path => $bam_path,
            model => $model,
            build => $build,
        };
    }

    return sort { $a->{name} cmp $b->{name} } @samples;
}


sub status {
    print STDERR "Status...\n";
    return _status ( _load_samples() );
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
