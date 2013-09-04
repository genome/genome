package Genome::Model::Tools::Predictor::Keggscan::BladeBlastBatcher;

use strict;
use warnings;

use Genome;

# The only place this module seems to be used is in Genome::Model::Tools::Predictor::Keggscan
# (which may not be used, either).  It used to be a complete cut-and-paste copy of
# Genome::Model::GenePrediction::Command::Pap::Blast::BladeBlastBatcher.  Until there's
# some way to know that this is not used anymore, we'll just pretend to be it.

class Genome::Model::Tools::Predictor::Keggscan::BladeBlastBatcher {
    is => 'Genome::Model::GenePrediction::Command::Pap::Blast::BladeBlastBatcher'
};

1;
