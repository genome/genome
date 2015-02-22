package Genome::Model::Tools::Sx::Trim::Flexbar::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Trim::Flexbar::Base {
    is => 'Genome::Model::Tools::Sx::ExternalCmdBase',
    is_abstract => 1,
};

sub cmd_display_name { 'flexbar '.$_[0]->version; }

sub execute {
    my $self = shift;
    $self->debug_message('Run '.$self->cmd_display_name.'...');

    my $cmd = $self->build_command;
    return if not $cmd;

    my $output = $self->_init_ouptut;
    return if not $output;

    my $rv = $self->_run_command($cmd);
    return if not $rv;
    $self->debug_message('Run flexbar...OK');

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
        $self->error_message('Failed to open reader for flexbar output!');
        return;
    }

    $self->debug_message('Processing flexbar output...');
    while ( my $seqs = $output_reader->read ) {
        $output->write($seqs);
    }

    $self->_rm_tmpdir;

    $self->debug_message('Run flexbar...OK');
    return 1;
}

sub _required_type_and_counts_for_inputs {
    return ( 'sanger', [qw/ 1 2 /], );
}

sub build_command {
    my $self = shift;

    my $cmd = $self->SUPER::build_command;
    return if not $cmd;

    my @input_params = $self->_resolve_input_params;
    return if not @input_params;


    my $input_param_name = $self->input_param_name;
    $cmd .= ' --'.$input_param_name.' '.$input_params[0]->{file};
    if ( $input_params[1] ) {
        $cmd .= ' --'.$input_param_name.'2 '.$input_params[1]->{file};
    }
    $cmd .= ' --format '.$self->input_format_param_value;
    $cmd .= ' --target '.$self->_tmpdir.'/output.fastq',

    return $cmd;
}

1;

