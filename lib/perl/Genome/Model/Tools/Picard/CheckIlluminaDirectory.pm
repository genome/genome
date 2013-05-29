
package Genome::Model::Tools::Picard::CheckIlluminaDirectory;

# http://picard.sourceforge.net/command-line-overview.shtml#CheckIlluminaDirectory

use strict;
use warnings FATAL => 'all';

use Genome;

class Genome::Model::Tools::Picard::CheckIlluminaDirectory {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        basecalls_directory => {
            is  => 'String',
            doc => 'The basecalls output directory.',
        },
        lane => {
            is  => 'String',
            doc => 'Lane number.',
        },
        read_structure => {
            is  => 'String',
            doc => 'A description of the logical structure of clusters in an Illumina Run',
        },
    ],
};

sub help_brief {
    'Check that the files in an Illumina run directory are available, exist, and are reasonably sized for every tile/cycle'
}

sub help_detail {
    return <<EOS
    Check an Illumina run directory.  For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CheckIlluminaDirectory
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path . '/CheckIlluminaDirectory.jar';
    unless (-e $jar_path) {
        die('Failed to find '. $jar_path .'!  This command may not be available in version '. $self->use_version);
    }

    my %map_args = qw(
        basecalls_directory basecalls_dir
        lane lanes
    );

    my $args = join(' ',
        map {
            my $value = $self->$_;
            my $arg = $map_args{$_} || $_;
            defined($value) ? (uc($arg) . "='$value'") : ()
        } sort qw(
            basecalls_directory
            lane
            read_structure
        )
    );

    my $cmd = $jar_path . " net.sf.picard.illumina.CheckIlluminaDirectory $args";
    $self->run_java_vm(
        cmd          => $cmd,
    );
    return 1;
}

1;
__END__

