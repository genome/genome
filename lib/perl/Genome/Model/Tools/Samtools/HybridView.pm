package Genome::Model::Tools::Samtools::HybridView;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
use POSIX;
our $VERSION = '0.01';
use Cwd;
use File::Basename;
use File::Path;

class Genome::Model::Tools::Samtools::HybridView {
    is => 'Command',
    has_optional_input => [
    ref_fasta => {
        is=>'Text',
        is_optional=>0,
    },
    output_glf => {
        is=>'Text',
        is_optional=>0,
        is_output=>1,
    },
    bam => {
        is=>'Text',
        is_optional=>0,
    },
    
    ],

};

sub help_brief {
    "simulates reads and outputs a name sorted bam suitable for import into Genome::Model"
}

sub help_detail {
}
#/gscuser/dlarson/src/samtools-0.1.7a-hybrid/samtools-hybrid view -uh /gscmnt/ams1182/info/model_data/2879948110/build114973139/alignments/114811720.bam | /gscuser/dlarson/src/samtools-0.1.7a-hybrid/samtools-hybrid calmd -Aur - /gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa 2>  /dev/null | /gscuser/dlarson/src/samtools-0.1.7a-hybrid/samtools-hybrid pileup - -g -f  /gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa>  H_ME-20000492_2-20000492_2_a.1.glf
sub execute {
    $DB::single=1;
    my $self=shift;
    my $samtools_cmd = "/gscuser/dlarson/src/samtools-0.1.7a-hybrid/samtools-hybrid";
    my $bam = $self->bam;
    my $ref = $self->ref_fasta;
    my $view_cmd = "$samtools_cmd view -uh $bam";
    my $calmd_cmd = "$samtools_cmd calmd -Aur - $ref 2> /dev/null";
    my $pileup_cmd = "$samtools_cmd pileup - -g -f $ref";
    my $output = $self->output_glf;
    my $cmd = "$view_cmd | $calmd_cmd | $pileup_cmd > $output";
    my $rv = Genome::Sys->shellcmd(cmd=> $cmd, input_files=>[$bam, $ref]);
    if($rv != 1) {
        return;
    }
    return 1;
}

1;
