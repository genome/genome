package Genome::Model::Tools::Bed::CoverageBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::CoverageBed {
    is => ['Command'],
    has_input => [
        bed_file => {
            is => 'Text',
            doc => 'The BED format file of intervals to pull out bases from reference.',
            is_input=>1,
        },
        bam_file => {
            is=>'Text',
            is_input=>1,
            is_optional=>0,
        },
        output_file=> {
            is=>'Text',
            is_output=>1,
            is_input=>1,
            is_optional=>0,
        },
    ],
};


sub execute {
    my $self = shift;
   
    #coverageBed -abam /gscmnt/gc7001/info/model_data/2880954593/build116975444/alignments/116991607.bam -b CleftPalate03110402capture_chip_liftover.bed -d  > test_coverageBed.output
    my $bam = $self->bam_file;
    my $bed = $self->bed_file;
    unless(-s $bam && -s $bed) {
        $self->error_message("you supplied a faulty file. you need to think about your life choices and how you got to this point.");
        return;
    }
    my $output_file = $self->output_file;
    my $cmd = "coverageBed -abam $bam -b $bed -d | bgzip -c > $output_file";
    Genome::Sys->shellcmd(cmd=>$cmd);


    return 1;
}
