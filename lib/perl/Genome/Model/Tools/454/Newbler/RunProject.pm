package Genome::Model::Tools::454::Newbler::RunProject;

use strict;
use warnings;

class Genome::Model::Tools::454::Newbler::RunProject {
    is => 'Genome::Model::Tools::454::Newbler',
    has => [
            dir => {
                    is => 'String',
                    doc => 'pathname of the output directory for project',
                },
        ],
    has_optional => [
                     params => {
                                 is => 'String',
                                 doc => 'command line params to pass to newbler',
                             },
                 ],

};


sub help_brief {
    "launch assembly (or reference alignment) on a given newbler project using the regions which have been added"
}

sub help_detail {
    return <<"EOS"

EOS
}

sub execute {
    my $self = shift;

    $DB::single = $DB::stopper;
    my $params = $self->params || '';
    my $cmd = $self->full_bin_path('runProject') .' '. $params .' '. $self->dir;
    $self->status_message("Running: $cmd");
    my $rv = system($cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero return status from command '$cmd'");
        return
    }
    return 1;
}

1;

