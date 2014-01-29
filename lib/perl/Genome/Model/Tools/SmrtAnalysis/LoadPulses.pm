package Genome::Model::Tools::SmrtAnalysis::LoadPulses;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::LoadPulses {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        cmp_hdf5_file => {
            doc => 'The aligned reads in the format of cmp.h5',
            is_output => 1,
        },
        input_fofn => {
            doc => 'The input fofn of pls.h5 of bas.h5 files.',
        },
        skip => {
            is => 'Text',
            doc => 'This flag(any value other than zero is considered true) is used to skip this step.  Mainly a hack for CCS.',
            is_optional => 1,
        },
    ],
    has_optional_input => {
        metrics => {
            is => 'Text',
            doc => 'A comma delimited list of metrics to transfer from the pls.h5 or bas.h5 to the cmp.h5:  QualityValue, ClassifierQV, StartTime, WidthInFrames, PulseWidth, pkmid, IPD, Light, PreBaseFrames, WhenStarted, StartTimeOffset',
        },
    },
};


sub help_brief {
    'Load pulse information and quality values into a Compare file.'
}


sub help_detail {
    return <<EOS
loadPulses - Load pulse information and quality values into a Compare file
usage: loadPulses movieFile cmpFile [-metrics m1,m2,...] 
movieFile may be a movie file or a fofn of movie file names.
metrics m1,m2,... is a comma-separated list of metrics to print to the pulse file.
valid metrics are: 
    QualityValue,ClassifierQV,StartTime, 
    WidthInFrames,PulseWidth,pkmid,IPD,
    Light, PreBaseFrames
    WhenStarted, StartTimeOffset
By default, QualityValue, ClassifierQV, StartTime, PulseWidth,
  WidthInFrames, pkmid, and IPD are added
Using hdf version 1.8.4
EOS
}

sub execute {
    my $self = shift;

    if ($self->skip) {
        $self->debug_message('Skipping LoadPulses!');
        return 1;
    }

    my $cmd = $self->analysis_bin .'/loadPulses '. $self->input_fofn .' '. $self->cmp_hdf5_file;
    if (defined($self->metrics)) {
        $cmd .= ' -metrics '. $self->metrics;
    }
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_fofn],
    );
    return 1;
}


1;
