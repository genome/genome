package Genome::Model::Tools::Tabix::Index;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Tabix::Index {
    is => 'Genome::Model::Tools::Tabix',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The file to index',
            shell_args_position => 1,
        },
    ],
    has_optional_input => [
        preset => {
            is => 'Text',
            doc => 'Use one of the preset file formats',
            valid_values => ['vcf', 'gff', 'bed', 'sam', 'vcf', 'psltbl'],
        },
        sequence_column => {
            is => 'Integer',
            doc => 'Sequence name column',
        },
        start_column => {
            is => 'Integer',
            doc => 'Start position column',
        },
        end_column => {
            is => 'Integer',
            doc => 'End position column',
        },
        skip_lines => {
            is => 'Integer',
            doc => 'Skip first N lines',
        },
        comment_char => {
            is => 'Text',
            doc => 'Symbol for comment/meta lines',
            default_value => '#',
        },
        force => {
            is => 'Boolean',
            doc => 'Force overwriting of previous index',
            default_value => 0,
        },
    ],
};

sub help_brief {
    return "Index one or more regions from a tabix indexed file.";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt tabix fetch in.vcf.gz --regions 1:180000-190000 --output-file=out.vcf --print-header
EOS
}

sub get_options {
    my $self = shift;
    my %pmap = (
        'preset' => 'p',
        'sequence_column' => 's',
        'start_column' => 'b',
        'end_column' => 'e',
        'skip_lines' => 'S',
        'comment_char' => 'c',
    );

    return join(" ",
            map {"-$pmap{$_}".$self->$_}
            grep {defined $self->$_} keys %pmap
        ) . ($self->force == 1 ? " -f" : "");
}

sub execute {
    my $self = shift;
    my $tabix = $self->tabix_path;
    my $input = $self->input_file;
    my $opts = $self->get_options;

    my $cmd = "$tabix $opts $input";
    $self->status_message("Creating tabix index for file $input");
    Genome::Sys->shellcmd(cmd => $cmd);

    return 1;
}

1;
