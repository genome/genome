package Genome::Model::Tools::454::Newbler::RunAssembly;

use strict;
use warnings;

class Genome::Model::Tools::454::Newbler::RunAssembly {
    is => 'Genome::Model::Tools::454::Newbler',
    has => [
        ],
    has_optional => [
                     params => {
                                is => 'String',
                                is_param => 1,
                                doc => 'parameters to pass to newbler',
                            },
                     sff_files => {
                                   is => 'String',
                                   is_input => 1,
                                   is_many => 1,
                                   doc => 'a comma separated list of sff file names to align',
                               },
                     sff_dir => {
                                 is => 'String',
                                 is_input => 1,
                                 doc => 'a directory of sff files to align',
                             },
                     assembly_dir => {
                                     is => 'String',
                                     is_input => 1,
                                     doc => 'the directory to perform the assembly/alignment',
                                 },
                 ],

};


sub help_brief {
    "launch a NON-PROJECT-BASED newbler denovo assembly"
}

sub help_detail {
    return <<"EOS"

EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    unless ($self->sff_files || $self->sff_dir) {
        $self->error_message("Must provide either command line option sff_files and/or sff_dir");
        return;
    }
    return $self;
}

sub execute {
    my $self = shift;

    $DB::single = $DB::stopper;

    my $params = $self->params || '';
    if ($self->assembly_dir) {
        $params .= ' -o '. $self->assembly_dir;
    }
    my $cmd = $self->full_bin_path('runAssembly') .' '. $params;
    if ($self->sff_dir) {
        $cmd .= ' '. $self ->sff_dir;
    }
    if ($self->sff_files) {
        for my $sff_file ($self->sff_files) {
            $cmd .= ' '. $sff_file;
        }
    }
    my $rv = system($cmd);
    $self->status_message('Running: '. $cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero return status from command '$cmd'");
        return
    }
    return 1;
}

1;

