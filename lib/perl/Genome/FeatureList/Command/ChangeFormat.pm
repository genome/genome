package Genome::FeatureList::Command::ChangeFormat;

use strict;
use warnings;

use Genome;

class Genome::FeatureList::Command::ChangeFormat{
    is => 'Command::V2',
    has_input => [
        feature_list => { is => 'Genome::FeatureList', doc => 'The feature list whose status is to be changed', shell_args_position => 1},
        format => { is => 'VARCHAR2', len => 64, doc => 'The format that the feature_list will be given', valid_values => ['1-based', 'true-BED', 'multi-tracked', 'multi-tracked 1-based'], },
    ],
};

sub help_brief {
   "Change the format of a feature-list.";     
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt feature-list change-format --feature-list 'example' --format 'true-BED'
EOS
}


sub help_detail {                           
    return <<EOS 
Updates the format of a feature-list.  This tool will currently only change the
status from 'unknown' to one of the other valid statuses.
EOS
}

sub execute {
    my $self = shift;
    my $feature_list = $self->feature_list;

    if ($feature_list->format ne 'unknown'){
        $self->error_message("feature-list " .  $feature_list->name . " does not have format set to 'unknown'");
        return;
    }

    $feature_list->format($self->format);
    return 1;
}

1;
