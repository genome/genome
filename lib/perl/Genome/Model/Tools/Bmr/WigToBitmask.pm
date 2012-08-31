package Genome::Model::Tools::Bmr::WigToBitmask;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Bit::Vector;

class Genome::Model::Tools::Bmr::WigToBitmask
{
  is => 'Command',
  has => [
  wig_file => {
    type => 'String',
    is_optional => 0,
    doc => 'The wiggle file to convert. Must be the fixedstep format with step=1 and span=1',
  },
  output_file => {
    type => 'String',
    is_optional => 0,
    doc => 'The resulting bitmask file to write to disk',
  },
  reference_index => {
    type => 'String',
    is_optional => 1,
    default => '/gscmnt/gc2106/info/medseq/ckandoth/refseq/all_sequences.fa.fai',
    doc => 'The Samtools index of the reference genome for which you are creating a bitmask',
  },
  _bitmask => {
    type => 'HashRef',
    is_optional => 1,
    doc => 'This is an internally used hashref representation of the bitmask',
  },
  ]
};

sub help_brief
{
  return 'This reads in a WIG format file and converts it to a bitmask representation of the whole genome.';
}

sub help_detail
{
  return <<HELP;
This script takes a file in WIG format, and uses it to set bits in a bit vector representation of
the genome. This may be useful for quickly querying whether locations or regions overlap with
others on a genome wide basis. By default, the bitmask is stored on the object. Thus calling this
script with no --output-file option from the command line will result in no output. Very sad. This
script was designed for coverage files, Thus non-zero values in the WIG file are flipped on in the
bitmask and 0 values are left off. If you want the reverse behavior, use the --exclude-regions
option. This script does not support track definitions nor the variableStep specification. See
https://gscweb.gsc.wustl.edu/wiki/File_Format/Wiggle_Track_Format_(WIG_-_WTF) for more information
on the wiggle format.
The bitmask itself is a hash reference to a hash with chromosome names as keys and the value a
Bit::Vector object where each position in the chromosome is represented by a bit. The length of
each chromosome is determined from the provided index file.
HELP
}

sub execute
{
  my $self = shift;
  #check inputs
  my $reference_index = $self->reference_index;
  my $wig_file = $self->wig_file;
  if((-z $reference_index) || (-z $wig_file))
  {
    $self->error_message("Input file missing or of zero size");
    return;
  }
  #read in reference sequences and initialize a bitvector hashed by chromosome
  my $ref_fh = IO::File->new($reference_index) or die "Couldn't open $reference_index. $!\n";
  my %genome = ();
  while(my $line = $ref_fh->getline)
  {
    my ($chr, $length) = split /\t/, $line;
    #Adding 1 to the length to allow indexing coordinates with a base of 1 instead of 0
    $genome{$chr} = Bit::Vector->new($length+1) or die "Failed to create bitvector. $!\n";
  }
  $ref_fh->close;

  my $wig_fh = IO::File->new($wig_file) or die "Couldn't open $wig_file. $!\n";
  my ($chromosome, $starting_origin, $block_idx) = (0, 0, 0);

  #Check the first line for Broad's "track" header to avoid having to check every line for it
  my $line = $wig_fh->getline;
  #If Broad's header wasn't used, then reset the file handle. Otherwise, proceed with next line
  if( $line !~ m/^track/ )
  {
    $wig_fh->close;
    $wig_fh = IO::File->new($wig_file);
  }
  while($line = $wig_fh->getline)
  {
    chomp $line;
    if($line =~ m/^fixed/)
    {
      ($chromosome, $starting_origin) = $line =~ /^fixedStep\s+chrom=(\S+)\s+start=(\d+)\s+/;
      $chromosome =~ s/chr//; #Remove a 'chr' prefix if any
      unless(exists($genome{$chromosome}))
      {
        $self->error_message("Chr$chromosome in line " . $wig_fh->input_line_number . " does not exist in reference file.");
        return;
      }
      $block_idx = $starting_origin;
    }
    else #Handle a line that we can safely assume contains a 0 or 1
    {
      if($line == 1)
      {
        $genome{$chromosome}->Bit_On($block_idx);
      }
      ++$block_idx; #We can also safely assume that step size is 1
    }
  }
  $wig_fh->close;
  $self->_bitmask(\%genome);
  if($self->output_file)
  {
    $self->write_genome_bitmask($self->output_file, $self->_bitmask);
  }
  return 1;
}

sub bitmask
{
  my $self = shift;
  return $self->_bitmask;
}

sub write_genome_bitmask
{
  my ($self, $filename, $genome_ref) = @_;
  unless($filename || $genome_ref)
  {
    $self->error_message("Missing arguments to write_genome_bitmask(). Cannot proceed.");
    return;
  }
  #do some stuff to write this to a file without making it suck
  my $out_fh = IO::File->new($filename,">:raw") or die "Unable to write to $filename. $!";
  my $header_string = join("\t", map {$_ => $genome_ref->{$_}->Size()} sort keys %$genome_ref);
  my $write_string = pack 'N/a*', $header_string;
  my $write_result = syswrite($out_fh, $write_string);
  unless(defined $write_result && $write_result == length($write_string))
  {
    $self->error_message("Error writing the header");
    return;
  }
  for my $chr (sort keys %$genome_ref)
  {
    #first write the length in bytes 
    my $chr_write_string = $genome_ref->{$chr}->Block_Read();
    $write_result = syswrite $out_fh, pack("N",length($chr_write_string));
    unless(defined $write_result || $write_result != 4)
    {
      $self->error_message("Error writing the length of chromosome $chr");
      return;
    }
    $write_result = syswrite $out_fh, $genome_ref->{$chr}->Block_Read();
    unless(defined $write_result || $write_result != length($chr_write_string))
    {
      $self->error_message("Error writing the header");
      return;
    }
  }
  $out_fh->close;
  return 1;
}

sub read_genome_bitmask
{
  my ($self, $filename) = @_;
  unless($filename)
  {
    $self->error_message("No filename provided.");
    return;
  }
  #do some stuff to read this from a file without making it suck
  my $in_fh = IO::File->new($filename,"<:raw")or die "Unable to read from $filename. $!";
  my $read_string;
  sysread $in_fh, $read_string, 4;
  my $header_length = unpack "N", $read_string;
  sysread $in_fh, $read_string, $header_length;
  my $header_string = unpack "a*",$read_string;
  my %genome = split /\t/, $header_string; #each key is the name, each value is the size in bits

  #now read in each one
  foreach my $chr (sort keys %genome)
  {
    $genome{$chr} = Bit::Vector->new($genome{$chr}) or die "Failed to create bitmask $!\n";
    sysread $in_fh, $read_string, 4;
    my $chr_byte_length = unpack "N", $read_string;
    my $chr_read_string;
    sysread $in_fh, $chr_read_string, $chr_byte_length;
    $genome{$chr}->Block_Store($chr_read_string);
  }
  $in_fh->close;
  $self->_bitmask(\%genome);
  return \%genome;
}

1;
