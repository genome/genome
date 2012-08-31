package Genome::Model::Tools::Music::Galaxy;

use strict;
use warnings;

use Genome;
use File::Basename;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Galaxy {
  is => "Genome::Model::Tools::Music::Play",
  has_input => [
    output_bundle => {
      is => 'Text',
      doc => 'Location where Galaxy would like the bundle of Music outputs to be saved',
    },
  ],
  has => [
    output_dir => {
      is_input => 0,
      is_output => 0,
      is_optional => 1,
      default => '',
    },
  ],
};

sub _is_hidden_in_docs { 1 }

sub execute {
  my $self = shift;

  my $output_dir = Genome::Sys->create_temp_directory();
  $self->output_dir($output_dir);

  $self->bam_list($self->create_new_bam_list);

  $self->SUPER::_execute_body(@_);

  my $tar_path = $self->output_bundle;
  my $cmd = "tar -cf $tar_path -C $output_dir .";
  my $rv = Genome::Sys->shellcmd(
      cmd => $cmd,
  );

  return 1;
}

sub create_new_bam_list {
  my $self = shift;
  my $original_bam_file = $self->bam_list;
  my $output_dir = Genome::Sys->create_temp_directory();
  my $new_bam_list = "$output_dir/bam_list";
  my @bams = Genome::Sys->read_file($original_bam_file);
  my $new_bam_fh = Genome::Sys->open_file_for_writing($new_bam_list);

  for (@bams) {
    chomp;
    my ($sample, $first_bam, $second_bam) = split("\t", $_);
    my $first_new_path = $self->index_and_link_bam_file($output_dir, $first_bam);
    my $second_new_path = $self->index_and_link_bam_file($output_dir, $second_bam);
    $new_bam_fh->print("$sample\t$first_new_path\t$second_new_path\n");
  }

  $new_bam_fh->close();

  return $new_bam_list;
}

sub index_and_link_bam_file {
  my $self = shift;
  my $output_dir = shift;
  my $bam = shift;

  my $filename = $output_dir . "/" . basename($bam);

  Genome::Sys->create_symlink($bam, $filename);
  Genome::Sys->shellcmd( cmd => "samtools index $filename");

  return $filename;
}

1;
