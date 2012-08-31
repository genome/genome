#:boberkfe this looks like a good place to use memcache to cache up some build status.
#:boberkfe when build events update, stuff their status into memcache.  gathering info otherwise
#:boberkfe can get reaaaaal slow.

package Genome::Model::Build::View::Status::Html;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::View::Status::Html {
    is => 'Genome::View::Status::Html',
};

1;
