use strict;
use warnings;
use Genome;

package Genome::Model::Command::Export::PublicationMethods;

class Genome::Model::Command::Export::PublicationMethods {
  is => 'Command::V2',
  has_many_input => [
    models => { is => 'Genome::Model', shell_args_position => 1,
            doc => 'the model or models to describe' },
  ],
  doc => 'output text for use in the methods section of a publication'
};

sub help_synopsis {
  return <<EOS

genome model export publication-methods 12345

genome model export publication-methods mymodelname

genome model export publication-methods "name like 'PNC2%'"

genome model export publication-methods "model_groups.name like 'hack%'"

EOS
}

sub help_detail {
  return <<EOS
This tool takes a model or group of models and outputs a publication description.

The model type in question must implement the publication_description method.

If the models which do implement that method return different descriptions, each
description will be shown, preceded by the associated list of models.
EOS
}

sub execute {
  my $self = shift;
  my @models = $self->models;
  my %descriptions;
  for my $model ($self->models) {
    unless ($model->can('publication_description')) {
      warn "model " . $model->__display_name__ . " does not implement the publication_description method!";
      next;
    }
    my $desc = $model->publication_description;
    my $list = $descriptions{$desc} ||= [];
    push @$list, $model;
  }
  my @d = keys %descriptions;

  if (@d == 1) {
    print "\n",$d[0],"\n\n";
    return 1;
  }
  elsif (@d == 0) {
    $self->error_message("No publication description(s) available.");
    return 0;
  }
  else {
    for my $d (sort @d) {
      my $model_list = $descriptions{$d};
      print "For:\n\t" . join("\n\t", map { $_->__display_name__ } @$model_list) . "\n\n";
      print "\n",$d,"\n";
    }
    print "\n";
    return 1;
  }
}

1;
