package Genome::Qc::Command::Config::Diff;

use strict;
use warnings 'FATAL';

use Genome;

use Text::Diff qw( diff );

class Genome::Qc::Command::Config::Diff {
    is => 'Genome::Command::Viewer',
    has => {
        qc_config_a => {
            is => 'Genome::Qc::Config',
            shell_args_position => 1,
            doc => 'QC config A to diff.',
        },
        qc_config_b => {
            is => 'Genome::Qc::Config',
            shell_args_position => 2,
            doc => 'QC config B to diff.',
        },
    },
    doc => 'Diff QC configurations',
};

sub help_detail { 'Diff QC configurations' }

sub write_report {
    my ($self, $width, $handle) = @_;
    my $stringa = $self->qc_config_a->config_to_yaml;
    my $stringb = $self->qc_config_b->config_to_yaml;
    printf ($handle "\nDiffing\n%s\n vs.\n%s\n", $self->qc_config_a->__display_name__, $self->qc_config_b->__display_name__);
    my $diff = diff \$stringa, \$stringb, { STYLE => "Context" };
    print $handle $diff;
}

1;

