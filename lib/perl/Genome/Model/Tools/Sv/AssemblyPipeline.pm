package Genome::Model::Tools::Sv::AssemblyPipeline;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Sv::AssemblyPipeline {
  is => 'Command',
  doc => 'Tools in the Assembly Pipeline for SV related output files.'
};

sub help_synopsis {
  my $self = shift;
  return <<HELP;
gmt sv assembly-pipeline ...
HELP
}

1;
