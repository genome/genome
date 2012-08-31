package Genome::Model::Tools::Cat;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Cat {
    is => 'Command',
    has => [
        dest => {
            doc => 'Destination file',
            shell_args_position => 1,
            is_input => 1,
            is_output => 1
        },
        source => {
            doc => 'List of source files',
            shell_args_position => '2',
            is_input => 1,
            is_many => 1
        },
    ]
};

sub help_brief {
    "Concatenate sources files into destination",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt cat dest source... 
EOS
}

sub help_detail {
    return <<EOS 
Concatenate is used in the Pindel workflow.
EOS
}

sub execute {
    my $self = shift;

    my $buf = '';
    my $n = 0;

    my $out = Genome::Sys->open_file_for_writing($self->dest);
    for my $file ($self->source) {
        my $in = Genome::Sys->open_file_for_reading($file);

        do {
            $n = $in->read($buf,4096,0);
            if (!defined $n) {
                $self->error_message('Failed to read: ' . $file);
                return;
            }

            $out->print($buf);
        } while ($n == 4096);

        $in->close;
    }

    $out->close;

    return 1;
}

1;
