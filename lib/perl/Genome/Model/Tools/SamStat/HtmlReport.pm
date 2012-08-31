package Genome::Model::Tools::SamStat::HtmlReport;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SamStat::HtmlReport {
    is => 'Genome::Model::Tools::SamStat::Base',
    has_input => [
        input_files => {},
    ],
    has_output => [
        output_files => {
            is_optional => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    my @input_files;
    if ( ref($self->input_files) eq 'ARRAY' ) {
        @input_files = @{$self->input_files};
    } else {
        @input_files = split(',', $self->input_files);
    }
    my $cmd = $self->samstat_path .' '. join(' ',@input_files);
    my @output_files = map { $_ .'.html' } @input_files;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => \@output_files,
    );
    $self->output_files(\@output_files);
    return 1;
}

