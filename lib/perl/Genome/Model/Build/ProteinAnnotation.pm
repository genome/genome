package Genome::Model::Build::ProteinAnnotation;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::ProteinAnnotation {
    is => 'Genome::Model::Build',
};

# method defined in the superclass for implementation.
sub files_ignored_by_diff {
    return qw(
        build.xml
        workflow.xml
        workflow.png
        keggscan/logs/\d+.*
        reports/Build_Initialized/report.xml
        reports/Build_Succeeded.*/report.xml
        interproscan/interpro.output.error
        keggscan/ancillary/genes.EC_only.+
    );
}

sub dirs_ignored_by_diff {
    return qw(
        logs/
        psortb/.*/
    );
}

sub create {
    die(__PACKAGE__ . ' is deprecated.');
}

1;
