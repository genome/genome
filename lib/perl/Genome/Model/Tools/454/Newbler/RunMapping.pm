package Genome::Model::Tools::454::Newbler::RunMapping;

use strict;
use warnings;

class Genome::Model::Tools::454::Newbler::RunMapping {
    is => 'Genome::Model::Tools::454::Newbler',
    has => [
            ref_seq => {
                        is => 'String',
                        is_input => 1,
                        doc => 'a reference genome to align with',
                    },
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
                     mapping_dir => {
                                     is => 'String',
                                     is_input => 1,
                                     doc => 'the directory to perform the mapping/alignment',
                                 },
                 ],

};


sub help_brief {
"gmt 454 newbler run-mapping --dir=DIR --ref-seq=PATH";
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
    unless (-s $self->ref_seq) {
        $self->error_message('The reference sequence provided does not exist or is zero size '. $self->ref_seq);
        return;
    }
    return $self;
}

sub execute {
    my $self = shift;

    $DB::single = $DB::stopper;

    my $params = $self->params || '';
    if ($self->mapping_dir) {
        $params .= ' -o '. $self->mapping_dir;
    }
    my $cmd = $self->full_bin_path('runMapping') .' '. $params .' '. $self->ref_seq;
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

