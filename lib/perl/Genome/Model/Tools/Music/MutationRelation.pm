package Genome::Model::Tools::Music::MutationRelation;

use warnings;
use strict;
use Carp;
use Genome;
use IO::File;
use POSIX qw( WIFEXITED );

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::MutationRelation {
  is => 'Command::V2',
  has_input => [
    bam_list => { is => 'Text', doc => "Tab delimited list of BAM files [sample_name, normal_bam, tumor_bam] (See Description)" },
    maf_file => { is => 'Text', doc => "List of mutations in MAF format" },
    mutation_matrix_file => { is => 'Text', doc => "Optionally store the sample-vs-gene matrix used during calculations.", is_optional => 1 },
    output_file => { is => 'Text', doc => "Results of mutation-relation tool", is_output => 1 },
    permutations => { is => 'Number', doc => "Number of permutations used to determine P-values", is_optional => 1, default => 100 },
    gene_list => { is => 'Text', doc => "List of genes to test, typically SMGs. If unspecified, all genes in MAF are tested.", is_optional => 1 },
    skip_non_coding => { is => 'Boolean', doc => "Skip non-coding mutations from the provided MAF file", is_optional => 1, default => 1 },
    skip_silent => { is => 'Boolean', doc => "Skip silent mutations from the provided MAF file", is_optional => 1, default => 1 },
  ],
  doc => "Identify relationships of mutation concurrency or mutual exclusivity in genes across cases.",
};

sub help_synopsis {
  return <<HELP
 ... music mutation-relation \\
        --maf-file /path/myMAF.tsv \\
        --permutations 1000 \\
        --output-file /path/mutation_relation.csv
HELP
}

sub help_detail {
  return <<HELP
This module parses a list of mutations in MAF format and attempts to determine relationships among
mutated genes. It employs a correlation test to see whether or not any two genes are mutated
concurrently (positive correlation) or mutually exclusively (negative correlation). Because of the
possibility of largely varying numbers of mutations present in different genes, P-values are
calculated using restricted permutations that take into account the distribution of mutation
counts among the samples. In the output file, 'pand' is the P-value for concurrent mutation
events, and 'pexc' is the P-value for mutually exclusive mutation events.
HELP
}

sub _additional_help_sections {
  return (
    "ARGUMENTS",
<<EOS

=over 4

=item --bam-list

=over 8

=item Provide a file containing sample names and normal/tumor BAM locations for each. Use the tab-
  delimited format [sample_name normal_bam tumor_bam] per line. This tool only needs sample_name,
  so all other columns can be skipped. The sample_name must be the same as the tumor sample names
  used in the MAF file (16th column, with the header Tumor_Sample_Barcode).

=back

=back

=cut

EOS
    );
}

sub _doc_authors {
  return <<EOS
 Nathan D. Dees, Ph.D.
 Qunyuan Zhang, Ph.D.
EOS
}

sub execute {
  # Parse input arguments
  my $self = shift;
  my $bam_list = $self->bam_list;
  my $maf_file = $self->maf_file;
  my $output_file = $self->output_file;
  my $gene_list = $self->gene_list;
  my $permutations = $self->permutations;
  my @all_sample_names; # names of all the samples, no matter if it's mutated or not
  my @genes_to_test; # the genes which will be tested for relationships
  my $skip_non_coding = $self->skip_non_coding;
  my $skip_silent = $self->skip_silent;

  # Parse out the names of the samples which should match the names in the MAF file
  my $sampleFh = IO::File->new( $bam_list ) or die "Couldn't open $bam_list. $!\n";
  my $line_count = 0;
  while( my $line = $sampleFh->getline )
  {
      $line_count++;
      next if ( $line =~ m/^#/ );
      chomp( $line );
      my ( $sample ) = split( /\t/, $line );
      if ($sample) {
          push( @all_sample_names, $sample );
      } else {
          warn("could not parse sample name from line $line_count of --bam_list");
      }
  }
  $sampleFh->close;

  # If user-specified, parse out the names of the genes that we are limiting our tests to
  if( defined $gene_list )
  {
      my $geneFh = IO::File->new( $gene_list ) or die "Couldn't open $gene_list. $!\n";
      while( my $line = $geneFh->getline )
      {
          next if ( $line =~ m/^#/ );
          chomp( $line );
          my ( $gene ) = split( /\t/, $line );
          push( @genes_to_test, $gene );
      }
      $geneFh->close;
  }

  # Create sample-gene matrix
  my $matrix_file = $self->create_sample_gene_matrix($maf_file, \@all_sample_names, \@genes_to_test, $skip_non_coding, $skip_silent);

  # Perform mutation-relation test using R
  my $R_cmd = "R --slave --args < " . __FILE__ . ".R $matrix_file $permutations $output_file";
  print "$R_cmd\n";
  WIFEXITED(system $R_cmd) or croak "Couldn't run: $R_cmd ($?)";

  return(1);
}

sub create_sample_gene_matrix {
    my $self = shift;
    my ( $maf_file, $all_sample_names_ref, $genes_to_test_ref, $skip_non_coding, $skip_silent ) = @_;
    my @all_sample_names = @{$all_sample_names_ref};
    my @genes_to_test = @{$genes_to_test_ref};

    # Create hash of mutations from the MAF file
    my %mutations;
    my %all_genes;

    # Parse the MAF file
    my $mafFh = IO::File->new( $maf_file ) or die "Couldn't open $maf_file. $!\n";
    while( my $line = $mafFh->getline )
    {
        next if( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $gene, $mutation_class, $sample ) = ( $cols[0], $cols[8], $cols[15] );

        # If the mutation classification is odd, quit with error
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
            next;
        }

        $all_genes{$gene}++;
        $mutations{$sample}{$gene}++;
    }
    $mafFh->close;

    # If the user specified a gene list, then check for genes that are not in the MAF
    if( scalar( @genes_to_test ) > 0 )
    {
        for( my $i = 0; $i < scalar( @genes_to_test ); ++$i )
        {
            unless( defined $all_genes{$genes_to_test[$i]} )
            {
                print "Skipping ", $genes_to_test[$i], " from specified gene-list since it was not found in the MAF file\n";
                splice( @genes_to_test, $i, 1 );
            }
        }
    }
    else
    {
        @genes_to_test = sort keys %all_genes;
    }

    # Write the input matrix to a file for use by the R code
    my $matrix_file;
    unless( $matrix_file = $self->mutation_matrix_file ) {
        $matrix_file = Genome::Sys->create_temp_file_path();
    }

    my $matrix_fh = new IO::File $matrix_file,"w";
    # Print input matrix file header
    my $header = join("\t","Sample",@genes_to_test);
    $matrix_fh->print("$header\n");

    # Print mutation relation input matrix
    for my $sample (sort @all_sample_names) {
        $matrix_fh->print($sample);
        for my $gene (@genes_to_test) {
            if (exists $mutations{$sample}{$gene}) {
                $matrix_fh->print("\t1");
            }
            else {
                $matrix_fh->print("\t0");
            }
        }
        $matrix_fh->print("\n");
    }

    return $matrix_file;
}

1;
