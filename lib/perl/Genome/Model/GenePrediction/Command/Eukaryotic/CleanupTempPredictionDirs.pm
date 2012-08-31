package Genome::Model::GenePrediction::Command::Eukaryotic::CleanupTempPredictionDirs;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::GenePrediction::Command::Eukaryotic::CleanupTempPredictionDirs {
    is => 'Command',
    has => [
        directory => {
            is => 'Path',
            is_input => 1,
            doc => 'All temp prediction directories located here are removed',
        },
    ],
};

sub help_brief { 
    return "Removes temp files and directories";
}

sub help_synopsis {
    return "Removes temp files and directories";
}

sub help_detail {
    return "Removes temp files and directories";
}

sub execute {
    my $self = shift;

    my @dirs = glob($self->directory . "/*_temp_predictions_*");
    $self->status_message("Found " . scalar @dirs . " temp directories to remove:\n" . join("\n", @dirs));

    for my $dir (@dirs) {
        $self->status_message("Removing $dir");
        next unless -d $dir;
        my $rv = system("rm -rf $dir");
        confess "Trouble removing $dir!" unless defined $rv and $rv == 0;
    }

    return 1;
}
     
1;

