package Genome::Info::UCSCConservation;

use strict;
use warnings;
use Genome;

my %ucsc_conservation_directories;

if ($ENV{GENOME_SYS_ID} and $ENV{GENOME_SYS_ID} ne 'GMS1') {
  %ucsc_conservation_directories =
  (
    36 => "/opt/gms/GMS1/fs/ams1100/info/v36_ucsc_conservation",
    37 => "/opt/gms/GMS1/fs/ams1100/info/v37_ucsc_conservation",
  );
}else{
  %ucsc_conservation_directories =
  (
    36 => "/gscmnt/ams1100/info/v36_ucsc_conservation",
    37 => "/gscmnt/ams1100/info/v37_ucsc_conservation",
  );
}

sub ucsc_conservation_directories{
    return %ucsc_conservation_directories;
}

1;
