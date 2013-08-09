package Genome::Model::Tools::SamStat::HtmlReport;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SamStat::HtmlReport {
    is => 'Genome::Model::Tools::SamStat::Base',
    has_input => [
        input_files => {},
        fix_xp => {
            is  => 'Boolean',
            doc => 'The flag to remove the XP tag (split alignment) in bam file.',
            default_value => 0,
            is_optional   => 1,
        },
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
    } 
    else {
        @input_files = split(',', $self->input_files);
    }

    my $samstat     = $self->samstat_path;
    my $input_files = join ' ', @input_files;
    my $cmd = $samstat .' '. $input_files;

    if ($self->fix_xp) {
        unless (@input_files == 1) {
            die $self->error_message("For now fix_xp option only applies to single bam input file");
        }
        $cmd = "samtools view -h $input_files " . '| sed \'s;\tXP:Z[^\t]*;;\' | '. "$samstat -f sam -n $input_files";
    }

    my @output_files = map { $_ .'.html' } @input_files;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files  => \@input_files,
        output_files => \@output_files,
    );

    $self->output_files(\@output_files);
    return 1;
}

