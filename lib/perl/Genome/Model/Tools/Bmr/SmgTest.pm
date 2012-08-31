package Genome::Model::Tools::Bmr::SmgTest;

use strict;
use warnings;
use IO::File;
use Genome;
use Cwd;

class Genome::Model::Tools::Bmr::SmgTest {
    is => 'Command',
    has => [
    gene_summary => {
        is => 'String',
        is_optional => 0,
        doc => 'File containing per-gene BMR info from GeneSummary.pm.',
    },
    output_file => {
        is => 'String',
        is_optional => 0,
        doc => 'The output file in which to store the pvalues from the SMG test.',
    },
    ]
};

sub help_brief {
    "Run the SMG test in R."
}

sub help_detail {
    "Takes as input output from gmt bmr gene-summary, and run's Qunyuan's SMG test which is coded in R."
}

sub execute {
    my $self = shift;
    my $rlibrary = "SMG_test.R";

    #Parse input
    my $genefile = $self->gene_summary;
    unless (-s $genefile) {
        $self->status_message("BMR file not found.\n");
        return;
    }

    #Call R for smg test
    my $test_output = $self->output_file;
    my $smg_test_cmd = "smg_test(in.file='$genefile',test.file='$test_output');";
    my $smg_test_rcall = Genome::Model::Tools::R::CallR->create(command=>$smg_test_cmd,library=>$rlibrary);
    $smg_test_rcall->execute;

    return 1;
}

1;
