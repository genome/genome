package Genome::Model::Tools::Sx::Trim::Far;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Trim::Far {
    is => 'Genome::Model::Tools::Sx::ExternalCmdBase',
    has => [
        __PACKAGE__->_cmd_properties,
        adapter => {
            is => 'Text',
            is_optional => 1,
            doc => 'Adaptor sequence to be removed.',
        },
        remove_revcomp => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Remove the reverse compl;ement of the single adapter.',
        },
        version => {
            is => 'Text',
            valid_values => [ __PACKAGE__->cmd_versions ],
            doc => 'Verion of FAR to use.',
        },
        _tmp_far_inputs => { is => 'Array', is_optional => 1, },
    ],
};

sub help_brief {
    return 'Trim reads with the [F]lexible [A]dapter [R]emover';
}

sub cmd_display_name {
    return 'FAR';
}

sub _cmd_versions {
    return (
        '1.7' => $ENV{GENOME_SW} . '/flexibleadapter/flexibleadapter-1.7/far',
        '1.84' => $ENV{GENOME_SW} . '/flexibleadapter/flexibleadapter-1.84/build/far',
        '2.0' => $ENV{GENOME_SW} . '/flexibleadapter/flexibleadapter-2.0/build/far',
        '2.17' => '/usr/bin/far2.17.0',
    );
}

sub cmd_versions {
    my $self = shift;
    my %versions = $self->_cmd_versions;
    return sort keys %versions;
}

sub _cmd_properties {
    return (
        adapters => {
            is => 'Text',
            is_optional => 1,
            doc => 'Fasta file of adapter sequences to be removed.',
        },
        adaptive_overlap => {
            is => 'Text',
            is_optional => 1,
            valid_values => [qw/ yes no /],
            default_value => 'no',
            doc => "If set to 'yes' the parameter '--min-overlap' will be treated as adaptive measure",
        },
        algorithm => {
            is => 'Text',
            is_optional => 1,
            valid_values => [qw/ needleman needlemanQuality /],
            doc => "Algorithm to use:\n needleman (complete global overlap alignment via needleman wunsch)\n needlemanQuality (quality based global alignment via needleman wunsch quality alignment)",
        },
        cut_off => {
            is => 'Text',
            is_optional => 1,
            default_value => 2,
            doc => 'Max number of allowed mismatches + indels for alignment per 10 bases.',
        },
        min_readlength => {
            is => 'Text',
            is_optional => 1,
            default_value => 18,
            doc => 'Minimum readlength in basepairs after adapter removal or read will be discarded.',
        },
        max_uncalled => {
            is => 'Text',
            is_optional => 1,
            default_value => 0,
            doc => 'Number of allowed uncalled bases in a read.',
        },
        min_overlap => {
            is => 'Text',
            is_optional => 1,
            default_value => 10,
            doc => 'minimum required overlap of adapter and sequence in basepairs.',
        },
        write_lengthdist => {
            is => 'Text',
            is_optional => 1,
            default_value => 'yes',
            doc => 'Writes a length distribution to the specified filename. By default it is: <targetfilename>.lengthdist (Options: yes/no) (default yes)',
        },
        nr_threads => {
            is => 'Text',
            is_optional => 1,
            default_value => 1,
            doc => 'Number of threads to use.',
        },
        trim_end => {
            is => 'Text',
            is_optional => 1,
            default_value => 'right',
            valid_values => [qw/ right left any left_tail right_tail /],
            doc => 'Decides on which end adapter removal is performed.',
        },
    );
}

sub execute {
    my $self = shift;

    my $resolve_adapters = $self->_resolve_adapters;
    return if not $resolve_adapters;

    my @input_params = $self->_resolve_input_params;
    return if not @input_params;

    my $output = $self->_init_ouptut;
    return if not $output;

    $self->debug_message('Run FAR...');
    my $cmd = $self->build_command;
    $cmd .= ' --source '.$input_params[0]->{file};
    $cmd .= ' --source2 '.$input_params[1]->{file} if $input_params[1];
    $cmd .= ' --format fastq-sanger';
    $cmd .= ' --target '.$self->_tmpdir.'/output.fastq',
    my $rv = $self->_run_command($cmd);
    return if not $rv;
    $self->debug_message('Run FAR...OK');

    my @fastq_files = glob $self->_tmpdir.'/output*.fastq';
    my @outputs = grep { $_ !~ /single/ } @fastq_files;
    if ( not @outputs ) {
        $self->error_message('Failed to find output files! Files in output directory: ');
        return;
    }
    $self->debug_message('Output: '.join(' ', @outputs));
    
    my $output_reader = Genome::Model::Tools::Sx::Reader->create(
        config => [ map { $_.':type=sanger' } @outputs ],
    );
    if ( not $output_reader ) {
        $self->error_message('Failed to open reader for FAR output!');
        return;
    }

    $self->debug_message('Processing FAR output...');
    while ( my $seqs = $output_reader->read ) {
        $output->write($seqs);
    }

    $self->_rm_tmpdir;

    return 1;
}

sub _resolve_adapters {
    my $self = shift;

    if ( $self->adapters and $self->adapter ) {
        $self->error_message('Cannot specify both adpaters file (adapters) and single adapter (adapter) params!');
        return;
    }

    if ( $self->adapters ) { # using adapters file
        $self->error_message('Cannot specify both adpaters file (adapters) and single adapter (adapter) params!') if $self->adapter;
        $self->error_message('Cannot specify to remove reverse complement adapter (remove_revcomp) and adpaters file (adapters) params!') if $self->remove_revcomp;
        return 1 if -s $self->adapters;
        $self->error_message('Adapters file does not exist! '.$self->adapters);
        return;
    }

    my $adapters_file = $self->_tmpdir.'/adapters.fasta';
    $self->adapters($adapters_file);
    my $writer = Genome::Model::Tools::Sx::PhredWriter->create(file => $adapters_file);
    if ( not $writer ) {
        $self->error_message('Failed to create fasta writer!');
        return;
    }

    my $adapter = uc($self->adapter);
    $writer->write({ id => 'Adapter', seq => $adapter, });
    if ( $self->remove_revcomp ) {
        my $revcomp_adapter = reverse $adapter;
        $revcomp_adapter =~ tr/TCGA/AGCT/;
        $writer->write({ id => 'Adapter-Rev-Complement', seq => uc( $revcomp_adapter ), });
    }

    return 1;
}

sub _required_type_and_counts_for_inputs {
    return ( 'sanger', [qw/ 1 2 /], );
}

1;

