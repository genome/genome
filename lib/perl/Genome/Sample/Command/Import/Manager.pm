package Genome::Sample::Command::Import::Manager;

use strict;
use warnings;

use Genome;

use Data::Dumper;
use Genome::Sample::Command::Import;
use IO::File;
use Switch;
use YAML;

class Genome::Sample::Command::Import::Manager {
    is => 'Command::V2',
    doc => 'Manage importing a group of samples including importing instrument data and creating and building models.',
    has => [
        working_directory => {
            doc => 'Directory to read and write.',
        },
    ],
    has_many_optional => [
        functions => {
            is => 'Text',
            valid_values => [qw/ create_samples /],
            doc => '',
        },
    ],
    has_calculated => [
        sample_csv_file => {
            calculate_from => 'working_directory',
            calculate => sub{ my $working_directory = shift; return $working_directory.'/samples.csv'; },
            doc => 'CSV file of samples and attributes. A column called "name" is required. The name should be dash (-) separated values of the nomenclature, indivdual id and sample id.',
        },
        config_file => {
            calculate_from => 'working_directory',
            calculate => sub{ my $working_directory = shift; return $working_directory.'/config.yaml'; },
        },
    ],
    has_optional_transient => [
        config => { is => 'Hash', },
        namespace => { is => 'Text', },
        samples => { is => 'Hash', },
    ],
};

sub execute {
    my $self = shift;

    my $samples = $self->_load_samples;
    return if not $samples;

    for my $function ( $self->functions ) {
        $self->status_message("Run '$function'...");
        $function = '_'.$function;
        my $rv = $self->$function;
        return if not $rv;
        $self->status_message("Run '$function'...OK");
    }
    $self->status($samples);

    return 1;
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    if ( not -d $self->working_directory ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ working_directory /],
            desc => 'Working directory does not exist or is not aq directory!',
        );
        return @errors;
    }

    my $config_error = $self->_load_config;
    if ( $config_error ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ config_file /],
            desc => $config_error,
        );
        return @errors;
    }

    my $sample_csv_error = $self->_load_sample_csv_file;
    if ( $sample_csv_error ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ sample_csv_file /],
            desc => $sample_csv_error,
        );
        return @errors;
    }

    return;
}

sub _load_config {
    my $self = shift;

    my $config_file = $self->config_file;
    if ( not -s $config_file ) {
        return 'Config file does not exist! '.$config_file;
    }

    my $config = YAML::LoadFile($config_file);
    if ( not $config ) {
        return 'Failed to load config file! '.$config_file;
    }
    $self->config($config);

    my $nomenclature = $config->{sample}->{nomenclature};
    if ( not $nomenclature ) {
        return 'No nomenclature in config file! '.$config_file;
    }

    my $namespace = Genome::Sample::Command::Import->namespace_for_nomenclature($nomenclature);
    if ( not $namespace ) {
        return 'Could not get namespace from miporter. Please ensure there is a config for the '.$nomenclature.' nomenclature.';
    }
    $self->namespace($namespace);

    return;
}

sub _load_sample_csv_file {
    my $self = shift;

    my $sample_csv_file = $self->sample_csv_file;
    if ( not -s $sample_csv_file ) {
        return 'Sample csv file does not exist! '.$sample_csv_file;
    }

    my $sample_csv_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $sample_csv_file,
        separator => ',',
    );
    if ( not $sample_csv_reader ) {
        return 'Failed to open sample csv! '.$sample_csv_file;
    }

    my %headers_not_found = ( name => 1, original_data_path => 1, );
    for my $header ( @{$sample_csv_reader->headers} ) {
        delete $headers_not_found{$header};
    }

    if ( %headers_not_found ) {
        return 'No '.join(' ', map { '"'.$_.'"' } keys %headers_not_found).' column in sample csv! '.$self->sample_csv_file;
    }

    my %samples;
    while ( my $sample = $sample_csv_reader->next ) {
        my $name = delete $sample->{name};
        $samples{$name} = {
            name => $name,
            original_data_path => delete $sample->{original_data_path},
            attributes => [],
            sample_attributes => [],
            individual_attributes => [],
        };
        for my $attr ( sort keys %$sample ) {
            my $attr_key = 'attributes';
            if ( $attr =~ /^(sample|individual)\./ ) {
                $attr_key = $1.'_attributes';
            }
            push @{$samples{$name}->{$attr_key}}, $attr."=\'".$sample->{$attr}."\'";
        }
    }
    $self->samples(\%samples);

    return;
}

sub _load_samples {
    my $self = shift;

    my $set_job_status_to_samples = $self->_set_job_status_to_samples;
    return if not $set_job_status_to_samples;

    my $samples = $self->samples;
    my %instrument_data = map { $_->sample->name, $_ } Genome::InstrumentData::Imported->get(
        original_data_path => [ map { $_->{original_data_path} } values %$samples ],
        '-hint' => [qw/ sample bam_path /],
    );

    my %genome_samples = map { $_->name, $_ } Genome::Sample->get(
        'name in' => [ keys %$samples ],
        '-hint' => [qw/ models instrument_data /],
    );

    for my $name ( keys %$samples ) {
        $samples->{$name}->{sample} = $genome_samples{$name};
        $samples->{$name}->{inst_data} = $instrument_data{$name};
        $samples->{$name}->{bam_path} = eval{ $samples->{inst_data}->bam_path };
        $samples->{$name}->{model} = eval{ ($samples->{inst_data}->models)[-1]; }; #FIXME! needs to get only model for this situation
        $samples->{$name}->{build} = eval{ $samples->{model}->latest_build };
    }

    $self->samples($samples);
    return 1;
}

sub _set_job_status_to_samples {
    my $self = shift;

    my $samples = $self->samples;
    Carp::confess('Need samples to load job status!') if not $samples;

    return 1 if not $self->config->{'job dispatch'};

    my $job_list_cmd = $self->config->{'job dispatch'}->{list}->{'command'};
    if ( not $job_list_cmd ) {
        $self->error_message('No job list "command" in config! '.YAML::Dump($self->config));
        return;
    }

    my $name_column = $self->config->{'job dispatch'}->{list}->{'name column'};
    if ( not $name_column ) {
        $self->error_message('No job list "name column" in config! '.YAML::Dump($self->config));
        return;
    }
    $name_column--;

    my $status_column = $self->config->{'job dispatch'}->{list}->{'status column'};
    if ( not $status_column ) {
        $self->error_message('No job list "status column" in config! '.YAML::Dump($self->config));
        return;
    }
    $status_column--;

    $job_list_cmd .= ' 2>/dev/null |';
    my $fh = IO::File->new($job_list_cmd);
    while ( my $line = $fh->getline ) {
        my @tokens = split(/\s+/, $line);
        my $name = $tokens[ $name_column ];
        next if not $samples->{$name};
        $samples->{$name}->{job_status} = lc $tokens[ $status_column ];
    }

    return 1;
}

sub _create_samples {
    my $self = shift;

    my $samples = $self->samples;
    Carp::confess('Need samples to create!') if not $samples;

    my $namespace = $self->namespace;
    my $importer_class_name = Genome::Sample::Command::Import->importer_class_name_for_namespace($namespace);
    my @importer_property_names = Genome::Sample::Command::Import->importer_property_names_for_namespace($namespace);
    my %genome_samples = map { $_->name => $_ } Genome::Sample->get(name => [ keys %$samples ]);
    for my $sample ( values %$samples ) {
        my $genome_sample = $genome_samples{ $sample->{name} };
        if ( not $genome_sample ) {
            my %import_params;
            $import_params{name} = $sample->{name};
            $import_params{sample_attributes} = $sample->{sample_attributes} if $sample->{sample_attributes};
            $import_params{individual_attributes} = $sample->{individual_attributes} if $sample->{individual_attributes};
            my $importer = $importer_class_name->create(%import_params);
            if ( not $importer ) {
                $self->error_message('Failed to create sample importer for sample!');
                return;
            }
            if ( not $importer->execute ) {
                $self->error_message('Failed to execute sample importer for sample!');
                return;
            }
            $genome_sample = $importer->_sample;
            if ( not $genome_sample ) {
                $self->error_message('Executed the importer successfully, but failed to create sample!');
                return;
            }
        }
        $sample->{sample} = $genome_sample;
    }

    return 1;
}

sub get_sample_status {
    my ($self, $sample) = @_;
    Carp::confess('No sample to set status!') if not $sample;
    return 'sample_needed' if not $sample->{sample};
    return 'import_'.$sample->{job_status} if $sample->{job_status};
    return 'import_needed' if not $sample->{inst_data};
    return 'import_failed' if not defined $sample->{bam_path} or not -s $sample->{bam_path};
    return 'model_needed' if $sample->{inst_data} and not $sample->{model};
    return 'build_requested' if $sample->{model}->build_requested;
    return 'build_needed' if not $sample->{build};
    return 'build_'.lc($sample->{build}->status);
}

sub set_sample_status {
    my ($self, $sample) = @_;
    Carp::confess('No sample to set status!') if not $sample;
    return $sample->{status} = $self->get_sample_status($sample);
}

sub _status {
    my $self = shift;

    my $samples = $self->samples;
    Carp::confess('No samples to display status!') if not $samples;

    my %totals;
    my $status;
    for my $name ( sort { $a cmp $b } keys %$samples ) {
        my $sample = $samples->{$name};
        $totals{total}++;
        $self->set_sample_status($sample);
        $totals{ $sample->{status} }++;
        $totals{build}++ if $sample->{status} =~ /^build/;
        $status .= sprintf("%-20s %10s\n", $sample->{name}, $sample->{status});
    }
    print "$status\nSummary:\n".join("\n", map { sprintf('%-16s %s', $_, $totals{$_}) } sort { $a cmp $b } keys %totals)."\n";
    return 1;
}

sub cleanup_failed {
    my $self = shift;
    my @samples = _load_samples();
    for my $sample ( @samples ) {
        $self->set_sample_status($sample);
        next if $sample->{status} ne 'import_failed' or ( defined $sample->{bam_path} and -s $sample->{bam_path} );
        print $sample->{name}.' '.$sample->{status}." REMOVE!\n";
        $sample->{inst_data}->delete;
    }
    UR::Context->commit;
    return 1;
}

sub commands {
    my $self = shift;
    print STDERR "Commands...\n";
    my @samples = _load_samples();
    for my $sample ( @samples ) {
        $self->set_sample_status($sample);
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

1;

