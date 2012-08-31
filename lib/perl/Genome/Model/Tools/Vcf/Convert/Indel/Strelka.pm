package Genome::Model::Tools::Vcf::Convert::Indel::Strelka;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::Vcf::Convert::Indel::Strelka {
    is => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from strelka output'
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from strelka indel output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the indels.
HELP
}

sub source {
    my $self = shift;
    return "Strelka";
}

#Override most of the print header code so that Strelka's header is concatenated to the TGI default header
sub print_header{
    my $self = shift;
    my $source           = $self->source;
    my $public_reference = $self->_get_public_ref;

    my $input_fh = $self->_input_fh;
    my $output_fh = $self->_output_fh;

    $output_fh->print("##reference=$public_reference" . "\n");
    $output_fh->print("##phasing=none" . "\n");
    $output_fh->print("##Original Strelka header follows:" . "\n");

    while(my $line = <$input_fh>) {
      chomp($line);
      if ($line =~ /^\#\#/){
        $output_fh->print($line, "\n");
      }elsif($line =~ /^\#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+NORMAL\s+TUMOR/){
        last;
      }else{
        die "Bad header: $line";
      }
    }

    my @header_columns = $self->_get_header_columns;

    #column header:
    $output_fh->print( "#" . join("\t", @header_columns) . "\n");
    return 1;
}


#Override the entire conversion process so that Strelka VCF lines are passed through unaltered
# Loop through each input line, parse it, and print it to output
sub convert_file {
  my $self = shift;
  my $input_fh = $self->_input_fh;

  #Skip comments and check header line and die if we find a problem
  while(my $line = <$input_fh>) {
    chomp($line);
    if ($line =~ /^\#\#/){
      next;
    }elsif($line =~ /^\#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+NORMAL\s+TUMOR/){
      last;
    }else{
      last;
    }
  }

  #Simply send the data line as is without any parsing
  while(my $line = $self->get_record($input_fh)) {
    chomp $line;
    $self->write_line($line);
  }

  return 1;
}

