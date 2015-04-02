package Genome::Model::Tools::Picard::CheckIlluminaDirectory;

use strict;
use warnings FATAL => 'all';

use Genome;

class Genome::Model::Tools::Picard::CheckIlluminaDirectory {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        basecalls_directory => {
            is  => 'String',
            doc => 'The basecalls output directory',
            picard_param_name => 'BASECALLS_DIR',
        },
        lane => {
            is  => 'String',
            doc => 'Lane number',
            picard_param_name => 'LANES',
        },
        read_structure => {
            is  => 'String',
            doc => 'A description of the logical structure of clusters in an Illumina Run',
            picard_param_name => 'READ_STRUCTURE',
        },
    ],
};

sub help_brief {
    'Check that the files in an Illumina run directory are available, exist, and are reasonably sized for every tile/cycle'
}

sub help_detail {
    return <<EOS
    Check an Illumina run directory.  For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#CheckIlluminaDirectory
EOS
}

sub _jar_name {
    return 'CheckIlluminaDirectory.jar';
}

sub _java_class {
    return qw(picard illumina CheckIlluminaDirectory);
}

1;
