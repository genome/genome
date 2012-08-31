package Genome::Model::Tools::SmrtAnalysis::CmpH5AppendConsensus;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::CmpH5AppendConsensus {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        cmp_hdf5_file => { is_output => 1, },
        cons_hdf5_files => {},
    ],
};

sub create {
    my $class = shift;
    my %params = @_;
    my $cons_hdf5_files = delete($params{cons_hdf5_files});
    my $self = $class->SUPER::create(%params);
    unless ($self) { return; }
    $self->cons_hdf5_files($cons_hdf5_files);
    return $self;
}

sub execute {
    my $self = shift;

    my $cmd = '/gsc/scripts/gsc/techd/cmpH5AppendConsensus.py';
    my @input_files;
    if (defined($self->cons_hdf5_files)) {
        my @cons_hdf5_files = @{$self->cons_hdf5_files};
        my $cons_hdf5_files = join(',',@cons_hdf5_files);
        $cmd .= ' --consH5Files='. $cons_hdf5_files;
        push @input_files, @cons_hdf5_files;
    }
    if (defined($self->cmp_hdf5_file)) {
        $cmd .= ' --cmpH5File='. $self->cmp_hdf5_file;
        push @input_files, $self->cmp_hdf5_file;
    }
    $self->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
