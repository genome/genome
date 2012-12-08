package Genome::Model::ClinSeq::Command::Test;

#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::ClinSeq::Command::Test {
    is => 'Command::V2',
    has_input => [
        builds => { 
              is => 'Genome::Model::Build::SomaticVariation',
              is_many => 1,
              shell_args_position => 1,
              require_user_verify => 0,
              doc => 'somatic variation build(s) to summarize SVs from',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
        _messages => {
            is => 'HASH',
            doc => 'hash of messages',
            default => {},
        },

    ],
    doc => 'summarize the SVs of somatic variation build',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq test --outdir=/tmp/create_mutation_diagram/  'id in [129973671,129708625]'

EOS
}

sub help_detail {
    return <<EOS
Test code ideas 

EOS
}

sub execute {
  my $self = shift;

  #$self->append_message("1", "Test message 1");
  #$self->calculate_something();
  #$self->append_message("3", "Test message 3");

  #my %messages = %{$self->_messages};
  #print Dumper %messages;


  my @mutation_diagram_builds;
  push(@mutation_diagram_builds, Genome::Model::Build->get(129973671));
  push(@mutation_diagram_builds, Genome::Model::Build->get(129708625));

  my $mutation_diagram_cmd = Genome::Model::ClinSeq::Command::CreateMutationDiagrams->create(builds=>\@mutation_diagram_builds, outdir=>$self->outdir, collapse_variants=>1, max_snvs_per_file=>250, max_indels_per_file=>250);
  my $r = $mutation_diagram_cmd->execute();


  return 1;
}


sub append_message {
    my $self = shift;
    my $key = shift || die;
    my $value = shift || die;

    my $messages = $self->_messages;
    $messages->{$key} = $value;
    $self->_messages($messages);

    return 1;
}


sub calculate_something {
  my $self = shift;

  my $x = 2+2;
  my $message = "Test message 2.  The answer is $x";
  $self->append_message("2", $message);

  return 1;
}


1;


