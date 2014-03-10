package Genome::Model::Tools::Music::Proximity;

use warnings;
use strict;
use IO::File;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Proximity {
  is => 'Command::V2',
  has_input => [
    maf_file => { is => 'Text', doc => "List of mutations using TCGA MAF specifications v2.3" },
    output_dir => { is => 'Text', doc => "Directory where output files will be written" },
    max_proximity => { is => 'Text', doc => "Maximum allowed AA distance between 2 mutations", is_optional => 1, default => 7 },
    skip_non_coding => { is => 'Boolean', doc => "Skip non-coding mutations from the provided MAF file", is_optional => 1, default => 1 },
    skip_silent => { is => 'Boolean', doc => "Skip silent mutations from the provided MAF file", is_optional => 1, default => 1 },
  ],
  has_optional_output => [
    output_file => {is => 'Text', doc => "TODO" },
  ],
  doc => "Perform a proximity analysis on a list of mutations."
};

sub help_detail {
  return <<HELP
This module first calculates the amino acid position of each mutation in the MAF file within its
respective transcript. Then, for each mutation, two values are calculated: 1) the number of other
mutations on the same transcript within the proximity limit set by the max-proximity input
parameter, and 2) the distance to the closest other mutation in this nearby set. Only mutations
which have another mutation within close proximity are reported in the output-file.

In addition to the standard version 2.3 MAF headers, there needs to be 3 columns appended. These
column headers in the MAF must have these names in the header in order for the tool to find them:
   transcript_name - the transcript name, such as NM_000028
 amino_acid_change - the amino acid change, such as p.R290H
        c_position - the nucleotide position changed, such as c.869

The output is generated with the folowing column headers:
Mutations_Within_Proximity, Nearest_Mutation, Gene, Transcript, Affected_Amino_Acid(s), Chr,
Start, Stop, Ref_Allele, Var_Allele, Sample
HELP
}

sub help_synopsis {
  return <<EOS
 ... music proximity \\
        --maf-file input_dir/myMAF.tsv \\
        --output-dir output_dir/ \\
        --max-proximity 15
EOS
}

sub _doc_authors {
  return <<EOS
 Nathan D. Dees, Ph.D.
 Dan Koboldt, M.S.
 Cyriac Kandoth, Ph.D.
EOS
}

sub execute {
  my $self = shift;
  my $maf_file = $self->maf_file;
  my $output_dir = $self->output_dir;
  my $max_proximity = $self->max_proximity;
  my $skip_non_coding = $self->skip_non_coding;
  my $skip_silent = $self->skip_silent;
  $output_dir =~ s/(\/)+$//; # Remove trailing forward slashes if any

  # Check on all the input data before starting work
  print STDERR "MAF file not found or is empty: $maf_file\n" unless( -s $maf_file );
  print STDERR "Output directory not found: $output_dir\n" unless( -e $output_dir );
  return undef unless( -s $maf_file && -e $output_dir );

  # Output of this script will be written to this location in the output directory
  my $out_file = "$output_dir/proximity_report";
  $self->output_file($out_file);

  # Parse the header row in the MAF file
  my $maf_fh = IO::File->new( $maf_file ) or die "Couldn't open $maf_file. $!";
  my $maf_header = $maf_fh->getline;
  $maf_header = $maf_fh->getline while( $maf_header =~ /^#/ ); # Skip commented lines
  chomp( $maf_header );
  unless( $maf_header =~ /^Hugo_Symbol/ )
  {
      print STDERR "Could not find column headers in $maf_file\n";
      return undef;
  }

  # Check whether the required additional MAF columns were included in the MAF
  unless( $maf_header =~ m/c_position/ and $maf_header =~ m/amino_acid_change/ and $maf_header =~ m/transcript_name/ )
  {
      print STDERR "Could not find required additional columns in $maf_file\n";
      return undef;
  }

  # Find the indexes of all the MAF columns
  my $idx = 0;
  my %col_idx = map {($_, $idx++)} split( /\t/, $maf_header );

  # A hash to store statuses, and a hash to store variants and AA positions
  my %status;
  my $status = \%status;
  my %aa_mutations;

  # Load relevant data from MAF into hash
  while( my $line = $maf_fh->getline )
  {
      chomp $line;
      my @cols = split( /\t/, $line );

      # Fetch data from the generic MAF columns
      my ( $gene, $chr, $start, $stop, $mutation_class, $mutation_type, $ref_allele, $var_allele, $var2, $sample ) =
      ( $cols[0], $cols[4], $cols[5], $cols[6], $cols[8], $cols[9], $cols[10], $cols[11], $cols[12], $cols[15] );
      $var_allele = $var2 if ( $var_allele eq $ref_allele ); # Different centers interpret the 2 variant columns differently

      # Fetch data from the required additional MAF columns
      my ( $c_position, $aa_change, $transcript ) = ( $cols[$col_idx{c_position}], $cols[$col_idx{amino_acid_change}], $cols[$col_idx{transcript_name}] );

      # Create a key to uniquely identify each variant
      my $variant_key = join( "\t", $gene, $chr, $start, $stop, $ref_allele, $var_allele, $sample );


      #check that the mutation class is acceptable
      if( $mutation_class !~ m/^(Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Silent|Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|De_novo_Start_InFrame|De_novo_Start_OutOfFrame)$/ )
      {
          print STDERR "Unrecognized Variant_Classification \"$mutation_class\" in MAF file for gene $gene\n";
          print STDERR "Please use TCGA MAF Specification v2.3.\n";
          return undef;
      }

      # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
      if(( $skip_non_coding && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/ ) ||
          ( $skip_silent && $mutation_class =~ m/^Silent$/ ))
      {
          print "Skipping $mutation_class mutation in gene $gene.\n";
          $status{synonymous_mutations_skipped}++;
          next;
      }

      # Determine amino acid position and load into hash
      my @mutated_aa_positions = ();
      @mutated_aa_positions = $self->get_amino_acid_pos( $variant_key, $mutation_class, $c_position, $aa_change, $status );

      # Record data in hash if mutated aa position found
      if( scalar( @mutated_aa_positions ) > 0 )
      {
          push( @{$aa_mutations{$transcript}{$variant_key}{mut_AAs}}, @mutated_aa_positions );
      }
  }
  $maf_fh->close;

  # Evaluate proximity of mutated amino acids for each transcript
  for my $transcript ( keys %aa_mutations )
  {
      # For each variant hitting that transcript
      for my $variant ( keys %{$aa_mutations{$transcript}} )
      {
          # Initialize the search
          my @affected_amino_acids = @{$aa_mutations{$transcript}{$variant}{mut_AAs}};
          my $mutations_within_proximity = 0; # Variable for summing # of mutations within proximity
          my $min_proximity = $max_proximity + 1; # Current minimum proximity

          # For each OTHER variant hitting the transcript
          for my $other_variant ( keys %{$aa_mutations{$transcript}} )
          {
              # Ignore the current mutation
              next if $variant eq $other_variant;
              # Get affected amino acids from OTHER variant
              my @other_affected_amino_acids = @{$aa_mutations{$transcript}{$other_variant}{mut_AAs}};

              # Compare distances between amino acids
              my $found_close_one = 0;
              for my $other_variant_aa ( @other_affected_amino_acids )
              {
                  for my $variant_aa ( @affected_amino_acids )
                  {
                      my $distance = abs($other_variant_aa - $variant_aa);
                      # If distance is within range
                      if( $distance <= $max_proximity )
                      {
                          $found_close_one++;
                          $min_proximity = $distance if $distance < $min_proximity;
                      }
                  }
              }
              # Note that this variant is within proximity if applicable
              $mutations_within_proximity++ if $found_close_one;
          }
          # Now, save results in hash if there are any
          if ($mutations_within_proximity)
          {
              $aa_mutations{$transcript}{$variant}{muts_within_range} = $mutations_within_proximity;
              $aa_mutations{$transcript}{$variant}{min_proximity} = $min_proximity;
          }
      }
  }

# Print results
  my $out_fh = IO::File->new( $out_file, ">" ) or die "Couldn't open $out_file. $!";
  $out_fh->print( "Mutations_Within_Proximity\tNearest_Mutation\tGene\tTranscript\tAffected_Amino_Acid(s)\tChr\tStart\tStop\tRef_Allele\tVar_Allele\tSample\n" );
  for my $transcript ( keys %aa_mutations )
  {
      for my $variant ( keys %{$aa_mutations{$transcript}} )
      {
          if( exists $aa_mutations{$transcript}{$variant}{muts_within_range} )
          {
              my ( $gene, $chr, $start, $stop, $ref_allele, $var_allele, $sample ) = split( /\t/, $variant );
              my $affected_amino_acids = join( ",", sort @{$aa_mutations{$transcript}{$variant}{mut_AAs}} );
              my $line = join( "\t", $aa_mutations{$transcript}{$variant}{muts_within_range},
                  $aa_mutations{$transcript}{$variant}{min_proximity}, $gene, $transcript,
                  $affected_amino_acids, $chr, $start, $stop, $ref_allele, $var_allele, $sample );
              $out_fh->print( "$line\n" );
          }
      }
  }
  $out_fh->close();
  return 1;
}

################################################################################

=head2 get_amino_acid_pos

This subroutine deducts the amino acid position within the transcript using the c_position and amino_acid_position columns in the MAF.

=cut

################################################################################

sub get_amino_acid_pos {
    # Parse arguments
    my $self = shift;
    my ( $variant_key, $mut_class, $c_position, $aa_change, $status ) = @_;

    # Initialize variables
    my $tx_start = my $tx_stop = 0;
    my $aa_position_start = my $aa_position_stop = 0;
    my $inferred_aa_start = my $inferred_aa_stop = 0;
    my $aa_pos = my $inferred_aa_pos = 0;

    # Amino acid position determination
    if( $aa_change && $aa_change ne "NULL" && substr( $aa_change, 0, 1 ) ne "e" )
    {
        $aa_pos = $aa_change;
        $aa_pos =~ s/[^0-9]//g;
    }

    # Parse out c_position if applicable ##
    if( $c_position && $c_position ne "NULL" )
    {
        # If multiple results, parse both ##
        if( $c_position =~ '_' && !( $mut_class =~ 'splice' ))
        {
            ($tx_start, $tx_stop) = split( /\_/, $c_position );
            $tx_start =~ s/[^0-9]//g;
            $tx_stop =~ s/[^0-9]//g;

            if( $tx_stop < $tx_start )
            {
                $inferred_aa_start = $tx_stop / 3;
                $inferred_aa_start = sprintf( "%d", $inferred_aa_start ) + 1 if( $tx_stop % 3 ) ;
                $inferred_aa_stop = $tx_start / 3;
                $inferred_aa_stop = sprintf( "%d", $inferred_aa_stop ) + 1 if( $tx_start % 3 );
            }
            else
            {
                $inferred_aa_start = $tx_start / 3;
                $inferred_aa_start = sprintf( "%d", $inferred_aa_start ) + 1 if( $tx_start % 3 );
                $inferred_aa_stop = $tx_stop / 3;
                $inferred_aa_stop = sprintf( "%d", $inferred_aa_stop ) + 1 if( $tx_stop % 3 );
            }
        }
        else
        {
            my ( $tx_pos ) = split( /[\+\-\_]/, $c_position );
            $tx_pos =~ s/[^0-9]//g;
            $tx_start = $tx_stop = $tx_pos;
            if($tx_pos)
            {
                $inferred_aa_pos = $tx_pos / 3;
                $inferred_aa_pos = sprintf( "%d", $inferred_aa_pos ) + 1 if( $tx_pos % 3 );
                $inferred_aa_start = $inferred_aa_stop = $inferred_aa_pos;
            }
        }
    }

    # If we inferred aa start stop, proceed with it ##
    if( $inferred_aa_start && $inferred_aa_stop )
    {
        $aa_position_start = $inferred_aa_start;
        $aa_position_stop = $inferred_aa_stop;
        $status->{aa_position_inferred}++;
    }
    # Otherwise if we inferred aa position ##
    elsif( $aa_pos )
    {
        $aa_position_start = $aa_pos;
        $aa_position_stop = $aa_pos;
        $status->{c_position_not_available}++;
    }
    # Otherwise we were unable to infer the info ##
    else
    {
        $status->{aa_position_not_found}++;
        $self->debug_message( "Amino acid position not found for variant: $variant_key" );
        return;
    }

    # Proceed if we have aa_position_start and stop ##
    my %mutated_aa_positions;
    if( $aa_position_start && $aa_position_stop )
    {
        for( my $this_aa_pos = $aa_position_start; $this_aa_pos <= $aa_position_stop; $this_aa_pos++ )
        {
            $mutated_aa_positions{$this_aa_pos}++;
        }
    }

    my @mutated_aa_positions = keys %mutated_aa_positions;
    return @mutated_aa_positions;
}

1;
