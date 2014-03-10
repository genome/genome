package Genome::Model::Tools::Sx::Mask::Dust;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Mask::Dust{
    is => 'Genome::Model::Tools::Sx::ExternalCmdBase',
};

sub help_brief {
    return 'Mask reads with dust [BLAST]';
}

sub cmd_display_name {
    return 'dust';
}

sub executable_path { $ENV{GENOME_SW} . '/wu-blast/blast2x64_2006-05-04/filter/dust' }

sub execute {
    my $self = shift;

    my @input_params = $self->_resolve_input_params;
    return if not @input_params;

    my $output = $self->_init_ouptut;
    return if not $output;

    $self->debug_message('Run dust on '.@input_params.' files');
    my $cmd = $self->build_command;
    my @output_configs;
    for my $input_params ( @input_params ) {
        my $output_file = $input_params->{file}.'.dusted';
        my $output_config = "file=$output_file";
        $output_config .= ":qual_file=".$input_params->{qual_file} if $input_params->{qual_file};
        $output_config .= ":type=phred";
        push @output_configs, $output_config;

        $self->debug_message('Dust:   '.$input_params->{file});
        $self->debug_message('Output: '.$output_file);

        $cmd .= ' '.$input_params->{file}.' > '.$output_file;
        my $rv = $self->_run_command($cmd);
        return if not $rv;

        if ( not -s $output_file ) {
            $self->error_message('Failed to create output file!');
            return;
        }
    }
    $self->debug_message('Run dust...OK');

    my $output_reader = Genome::Model::Tools::Sx::Reader->create(config => \@output_configs);
    if ( not $output_reader ) {
        $self->error_message('Failed to open reader for dust output!');
        return;
    }

    $self->debug_message('Processing dust output...');
    while ( my $seqs = $output_reader->read ) {
        $output->write($seqs);
    }

    $self->_rm_tmpdir;

    return 1;
}

sub _required_type_and_counts_for_inputs {
    return ( 'phred', [ ], );
}

1;

