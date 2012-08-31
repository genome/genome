#!/usr/bin/env genome-perl

use strict;
use warnings;
use above 'Genome';

my $sample_name = $ARGV[0];

unless($sample_name) {
    print STDERR "Usage: dump_attribs.pl <sample_name>\n";
    exit(0);
}
my @samples;
if(my $temp = Genome::Sample->get(name=>"$sample_name")) {
  $samples[0]=$temp;
}

unless(scalar(@samples)==1) {
  my $mg = Genome::ModelGroup->get($ARGV[0]);
  unless($mg) {
    $mg = Genome::ModelGroup->get(name=>$ARGV[0]);
  }
  unless($mg) {
    print STDERR "Can't find $sample_name as a sample name, model group id, or model group name";
    exit(0);
  }
  my @models= $mg->models;
  @samples = sort {$a->name cmp $b->name } map {$_->subject} @models;
}

unless(@samples) {
    print STDERR "Sample or modelgroup not found!";
    exit(0);
}

for my $sample (@samples) {
    my @attribs = $sample->attributes;
    my $name = $sample->name;
    print "Sample: $name";
    for my $attrib (@attribs) {
        my $name = $attrib->attribute_label;
        my $value = $attrib->attribute_value;
        print "\t$name: $value";
    }
    print "\n";
}
