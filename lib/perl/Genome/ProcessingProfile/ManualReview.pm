package Genome::ProcessingProfile::ManualReview;

#:eclark 11/16/2009 Code review.

# Seems like there should better way to define a build that only does one thing.  Aside from that, this module has no problems (because it does nothing.)

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::ManualReview {
    is => 'Genome::ProcessingProfile::Staged',
};

sub stages {
    return (qw/
            manual_review
            /);
}

sub manual_review_job_classes {
    return (qw/
            Genome::Model::Event::Build::ManualReview::Run
        /);
}

sub manual_review_objects {
    return 1;
}

1;
