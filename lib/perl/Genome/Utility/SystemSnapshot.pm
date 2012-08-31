package Genome::Utility::SystemSnapshot;

use Genome; 

class Genome::Utility::SystemSnapshot {
    is  => 'Command',
};

sub new {
  my $class = shift;
  my $build_id = shift;
  my $datadir = shift;

  my $self = {
    build_id => $build_id,
    datadir  => $datadir,
    output   => "$datadir/output",
    errors   => "$datadir/errors",
    metrics  => "$datadir/$build_id.metrics.txt",
  };
  bless $self, $class;
  return $self;
}

sub pre_run {
  return;
}

sub run {
  return;
}

sub post_run {
  return;
}

sub report {
  my $self = shift;
  my $metrics;

  open S, "<$self->{metrics}" or die "Unable to open metrics file: $self->{metrics}: $!";
  my @lines = <S>;
  close S;
  foreach my $line (@lines) {
    chomp $line;
    next if ($line =~ /^$/);
    my ($metric,$value) = split(': ',$line);
    $metrics->{$metric} = $value;
  }
  return $metrics;
}

1;
