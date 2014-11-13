package Genome::Model::Build::ClinSeq::InputBuilds;

use Genome;
use strict;
use warnings;

class Genome::Model::Build::ClinSeq::InputBuilds {
    is => 'UR::Object',
};

sub resolve_somatic_builds{
  my $clinseq_build = shift;
  my $somatic_builds = shift;
  my $wgs_build = $clinseq_build->wgs_build;
  $somatic_builds->{$wgs_build->id}{build} = $wgs_build if $wgs_build;
  $somatic_builds->{$wgs_build->id}{type} = 'wgs' if $wgs_build;
  my $exome_build = $clinseq_build->exome_build;
  $somatic_builds->{$exome_build->id}{build} = $exome_build if $exome_build;
  $somatic_builds->{$exome_build->id}{type} = 'exome' if $exome_build;
}

sub resolve_rnaseq_builds {
  my $clinseq_build = shift;
  my $rnaseq_builds = shift;
  my $tumor_build = $clinseq_build->tumor_rnaseq_build;
  $rnaseq_builds->{$tumor_build->id}{build} = $tumor_build if $tumor_build;
  $rnaseq_builds->{$tumor_build->id}{type} = 'rnaseq' if $tumor_build;
  my $normal_build = $clinseq_build->normal_rnaseq_build;
  $rnaseq_builds->{$normal_build->id}{build} = $normal_build if $normal_build;
  $rnaseq_builds->{$normal_build->id}{type} = 'rnaseq' if $normal_build;
}

1;
