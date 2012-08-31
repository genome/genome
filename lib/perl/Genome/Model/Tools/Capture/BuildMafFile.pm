package Genome::Model::Tools::Capture::BuildMafFile;

##################################################################################################
# BuildMafFile - Generate MAF File using Capture Somatic Results
#
#  AUTHOR:    Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#  CREATED:  06/16/2010 by D.K.
#  MODIFIED:  08/13/2010 by D.K.
#
#  NOTES:
#
##################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome; # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##
my %stats = ();

class Genome::Model::Tools::Capture::BuildMafFile {
  is => 'Command',

  has => [                                # specify the command's single-value properties (parameters) <---
    data_dir    => { is => 'Text', doc => "Somatic-capture model build dir" , is_optional => 0},
    tumor_sample  => { is => 'Text', doc => "Name of the tumor sample" , is_optional => 0},
    normal_sample  => { is => 'Text', doc => "Name of the matched normal control" , is_optional => 0},
    output_file  => { is => 'Text', doc => "Output file to contain the MAF" , is_optional => 0},
    limit_variants  => { is => 'Text', doc => "List of variant positions to include" , is_optional => 1},
    build     => { is => 'Text', doc => "Reference genome build" , is_optional => 1, default => "36"},
    center     => { is => 'Text', doc => "Genome center name" , is_optional => 1, default => "genome.wustl.edu"},
    sequence_phase  => { is => 'Text', doc => "Sequencing phase" , is_optional => 1, default => "4"},
    sequence_source  => { is => 'Text', doc => "Sequence source" , is_optional => 1, default => "Capture"},
    sequencer  => { is => 'Text', doc => "Sequencing platform name" , is_optional => 1, default => "IlluminaGAIIx"},
  ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Build TCGA-friendly MAF file using somatic capture pipeline output"
}

sub help_synopsis {
    return <<EOS
Build TCGA-friendly MAF file using somatic capture pipeline output
EXAMPLE:  gmt capture build-maf-file
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS

EOS
}

################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
  my $self = shift;

  my $data_dir = $self->data_dir;
  my $tumor_sample = $self->tumor_sample;
  my $normal_sample = $self->normal_sample;
  my $output_file = $self->output_file;

  my $center = $self->center;
  my $build = $self->build;
  my $sequence_phase = $self->sequence_phase;
  my $sequence_source = $self->sequence_source;
  my $sequencer = $self->sequencer;


  ## Fix sample names ##
  my ($sample_name_header) = split(/\-/, $tumor_sample);
  
  ## If we have WashU sample names, correc tthem
  if(substr($sample_name_header, 0, 1) eq "H")
  {
    $normal_sample =~ s/$sample_name_header/TCGA/;
    $tumor_sample =~ s/$sample_name_header/TCGA/;    
  }

  ## Keep stats in a single hash ##

  my %stats = ();
  my %limit_variants = load_mutations($self->limit_variants) if($self->limit_variants);

  ## Check for required files in the data directory ##

  my $snv_tier1_file = $data_dir . "/" . "merged.somatic.snp.filter.novel.tier1";
  my $snv_annotation_file = $data_dir . "/" . "annotation.somatic.snp.transcript";
  my $indel_tier1_file = $data_dir . "/" . "merged.somatic.indel.filter.tier1";
  my $indel_annotation_file = $data_dir . "/" . "annotation.somatic.indel.transcript";
  my $gatk_indel_tier1_file = $data_dir . "/" . "gatk.output.indel.formatted.Somatic.tier1";
  my $gatk_indel_annotation_file = $data_dir . "/" . "annotation.somatic.gatk-indel.transcript";

  if(-e $snv_tier1_file && -e $snv_annotation_file)# && -e $indel_tier1_file && -e $indel_annotation_file)
  {
    ## Build or access the dbSNP file ##

    my $dbsnp_file = $snv_tier1_file . ".dbsnp";
    if(!(-e $dbsnp_file))
    {
      my $cmd = "gmt annotate lookup-variants --variant-file $snv_tier1_file --report-mode known-only --append-rs-id --output-file $snv_tier1_file.dbsnp";
      system($cmd);
    }

    ## Load dbSNPs ##

    my %dbsnp_rs_ids = load_dbsnps($snv_tier1_file . ".dbsnp");

    ## Open the outfile ##

    open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
    print OUTFILE join("\t", "Hugo_Symbol","Entrez_Gene_Id","GSC_Center","NCBI_Build","Chromosome","Start_position","End_position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","dbSNP_RS","dbSNP_Val_Status","Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode","Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2","Tumor_Validation_Allele1","Tumor_Validation_Allele2","Match_Norm_Validation_Allele1","Match_Norm_Validation_Allele2","Verification_Status","Validation_Status","Mutation_Status","Validation_Method","Sequencing_Phase","Sequence_Source","Score","BAM_file","Sequencer","chromosome_name_WU","start_WU","stop_WU","reference_WU","variant_WU","type_WU","gene_name_WU","transcript_name_WU","transcript_species_WU","transcript_source_WU","transcript_version_WU","strand_WU","transcript_status_WU","trv_type_WU","c_position_WU","amino_acid_change_WU","ucsc_cons_WU","domain_WU","all_domains_WU","deletion_substructures_WU") . "\n";
#    print OUTFILE join("\t", "Hugo_Symbol","Entrez_Gene_Id","GSC_Center","NCBI_Build","Chromosome","Start_position","End_position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","dbSNP_RS","dbSNP_Val_Status","Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode","Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2","Tumor_Validation_Allele1","Tumor_Validation_Allele2","Match_Norm_Validation_Allele1","Match_Norm_Validation_Allele2","Verification_Status","Validation_Status","Mutation_Status","Sequencing_Phase","Sequence_Source","Validation_Method","Score","BAM_file","Sequencer","chromosome_name_WU","start_WU","stop_WU","reference_WU","variant_WU","type_WU","gene_name_WU","transcript_name_WU","transcript_species_WU","transcript_source_WU","transcript_version_WU","strand_WU","transcript_status_WU","trv_type_WU","c_position_WU","amino_acid_change_WU","ucsc_cons_WU","domain_WU","all_domains_WU","deletion_substructures_WU") . "\n";

    ## Load the SNVs ##

    my %snvs = load_mutations($snv_tier1_file);
    my %snv_annotation = load_annotation($snv_annotation_file);

    foreach my $key (sort byChrPos keys %snvs)
    {
      $stats{'tier1_snvs'}++;

      my ($chromosome, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $key);
      my $snv = $snvs{$key};
      my $strand = "+";

      my @temp = split("\t", $snv);

      my $tumor_cns = $temp[4];
      $tumor_cns = $temp[12] if($temp[11] =~ '%');

      my $tumor_genotype = iupac_to_genotype($ref, $tumor_cns);
      my ($tumor_gt_allele1, $tumor_gt_allele2) = split(//, $tumor_genotype);

      if(!$self->limit_variants || $limit_variants{$key})
      {
        if($snv_annotation{$key})
        {
          $stats{'tier1_snvs_with_annotation'}++;
          my @annotation = split(/\t/, $snv_annotation{$key});
          my $gene_name = $annotation[6];

          ## Get the gene ID ##
          my $gene_id = 0;

          my @ea = Genome::Site::TGI::EntityAlias->get(alias => "$gene_name", alias_source => "HUGO", entity_type_name => "gene sequence tag");

          if(@ea)
          {
            my @tags = Genome::Site::TGI::SequenceTag->get(stag_id => [ map {$_->entity_id} @ea ]);
            if(@tags)
            {
              $gene_id = $tags[0]->ref_id;
            }
          }


          my $trv_type = $annotation[13];

          my $mutation_type = trv_to_mutation_type($trv_type);

          ## Get dbSNP Status
          my $dbsnp_rs = "novel";
          my $dbsnp_status = "unknown";

          if($dbsnp_rs_ids{$key})
          {
            $dbsnp_rs = $dbsnp_rs_ids{$key};
            $dbsnp_status = "unknown";
          }

          print OUTFILE join("\t", $gene_name,$gene_id,$center,$build,$chromosome,$chr_start,$chr_stop,$strand,$mutation_type,"SNP",$ref,$tumor_gt_allele1,$tumor_gt_allele2,$dbsnp_rs,$dbsnp_status,$tumor_sample,$normal_sample,$ref,$ref,"","","","","Unknown","Unknown","Somatic",$sequence_phase,$sequence_source,"","1","dbGAP",$sequencer, @annotation) . "\n";
          $stats{'tier1_snvs_written'}++;

        }
        else
        {
          warn "No annotation for $key in $snv_annotation_file!\n";
        }
      }
      else
      {
        $stats{'tier1_snvs_not_included'}++;
      }
    }


    ## Write indels to file ##
    my %indels_written = ();

    if(-e $indel_tier1_file && -e $indel_annotation_file)# && -e $indel_tier1_file && -e $indel_annotation_file)
    {
      my %indels = load_mutations($indel_tier1_file);
      my %indel_annotation = load_annotation($indel_annotation_file);

      foreach my $key (sort byChrPos keys %indels)
      {
        $stats{'tier1_indels'}++;

        my ($chromosome, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $key);
        my $indel = $indels{$key};
        my $strand = "+";

        my $variant_type = "Unknown";

        if($ref eq "0" || $ref eq "-" || length($var) > 1)
        {
          $variant_type = "INS";
        }
        else
        {
          $variant_type = "DEL";
        }

        my @temp = split("\t", $indel);

#        my $tumor_cns = $temp[4];
#        $tumor_cns = $temp[12] if($temp[11] =~ '%');

#        my $tumor_genotype = iupac_to_genotype($ref, $tumor_cns);
#        my ($tumor_gt_allele1, $tumor_gt_allele2) = split(//, $tumor_genotype);

        if(!$self->limit_variants || $limit_variants{$key})
        {
          if($indel_annotation{$key})
          {
            $stats{'tier1_indels_with_annotation'}++;
            my @annotation = split(/\t/, $indel_annotation{$key});

            my $tumor_gt_allele1 = $annotation[3];
            my $tumor_gt_allele2 = $annotation[4];

            my $gene_name = $annotation[6];

            ## Get the gene ID ##
            my $gene_id = 0;

            my @ea = Genome::Site::TGI::EntityAlias->get(alias => "$gene_name", alias_source => "HUGO", entity_type_name => "gene sequence tag");

            if(@ea)
            {
              my @tags = Genome::Site::TGI::SequenceTag->get(stag_id => [ map {$_->entity_id} @ea ]);
              if(@tags)
              {
                $gene_id = $tags[0]->ref_id;
              }
            }


            my $trv_type = $annotation[13];

            my $mutation_type = trv_to_mutation_type($trv_type);

            ## Get dbSNP Status
            my $dbsnp_rs = "novel";
            my $dbsnp_status = "unknown";

            if($dbsnp_rs_ids{$key})
            {
              $dbsnp_rs = $dbsnp_rs_ids{$key};
              $dbsnp_status = "unknown";
            }

            my $indel_key = "$chromosome\t$chr_start\t$chr_stop\t$variant_type";
            $indels_written{$indel_key} = 1;

            print OUTFILE join("\t", $gene_name,$gene_id,$center,$build,$chromosome,$chr_start,$chr_stop,$strand,$mutation_type,$variant_type,$ref,$tumor_gt_allele1,$tumor_gt_allele2,$dbsnp_rs,$dbsnp_status,$tumor_sample,$normal_sample,$ref,$ref,"","","","","Unknown","Unknown","Somatic",$sequence_phase,$sequence_source,"","1","dbGAP",$sequencer, @annotation) . "\n";
            $stats{'tier1_indels_written'}++;
          }
          else
          {
            warn "No annotation for $key in $indel_annotation_file!\n";
          }
        }
        else
        {
          $stats{'tier1_indels_not_included'}++;
        }
      }
    }


    ## Write GATK Indels to File ##

    if(-e $gatk_indel_tier1_file && -e $gatk_indel_annotation_file)# && -e $indel_tier1_file && -e $indel_annotation_file)
    {
      my %indels = load_mutations($gatk_indel_tier1_file);
      my %indel_annotation = load_annotation($gatk_indel_annotation_file);

      foreach my $key (sort byChrPos keys %indels)
      {
        $stats{'tier1_indels'}++;

        my ($chromosome, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $key);
        my $indel = $indels{$key};
        my $strand = "+";

        my $variant_type = "Unknown";

        if($ref eq "0" || $ref eq "-" || length($var) > 1)
        {
          $variant_type = "INS";
        }
        else
        {
          $variant_type = "DEL";
        }

        my @temp = split("\t", $indel);

#        my $tumor_cns = $temp[4];
#        $tumor_cns = $temp[12] if($temp[11] =~ '%');
#        my $tumor_genotype = iupac_to_genotype($ref, $tumor_cns);
#        my ($tumor_gt_allele1, $tumor_gt_allele2) = split(//, $tumor_genotype);

        if(!$self->limit_variants || $limit_variants{$key})
        {
          if($indel_annotation{$key})
          {
            $stats{'tier1_indels_with_annotation'}++;
            my @annotation = split(/\t/, $indel_annotation{$key});

            my $tumor_gt_allele1 = $annotation[3];
            my $tumor_gt_allele2 = $annotation[4];

            my $gene_name = $annotation[6];

            ## Get the gene ID ##
            my $gene_id = 0;

            my @ea = Genome::Site::TGI::EntityAlias->get(alias => "$gene_name", alias_source => "HUGO", entity_type_name => "gene sequence tag");

            if(@ea)
            {
              my @tags = Genome::Site::TGI::SequenceTag->get(stag_id => [ map {$_->entity_id} @ea ]);
              if(@tags)
              {
                $gene_id = $tags[0]->ref_id;
              }
            }


            my $trv_type = $annotation[13];

            my $mutation_type = trv_to_mutation_type($trv_type);

            ## Get dbSNP Status
            my $dbsnp_rs = "novel";
            my $dbsnp_status = "unknown";

            if($dbsnp_rs_ids{$key})
            {
              $dbsnp_rs = $dbsnp_rs_ids{$key};
              $dbsnp_status = "unknown";
            }

            my $indel_key = "$chromosome\t$chr_start\t$chr_stop\t$variant_type";

            if(!$indels_written{$indel_key})
            {
              $indels_written{$indel_key} = 1;
              print OUTFILE join("\t", $gene_name,$gene_id,$center,$build,$chromosome,$chr_start,$chr_stop,$strand,$mutation_type,$variant_type,$ref,$tumor_gt_allele1,$tumor_gt_allele2,$dbsnp_rs,$dbsnp_status,$tumor_sample,$normal_sample,$ref,$ref,"","","","","Unknown","Unknown","Somatic",$sequence_phase,$sequence_source,"","1","dbGAP",$sequencer, @annotation) . "\n";
              $stats{'tier1_indels_written'}++;
            }
          }
          else
          {
            warn "No annotation for $key in $indel_annotation_file!\n";
          }
        }
        else
        {
          $stats{'tier1_indels_not_included'}++;
        }
      }
    }





    close(OUTFILE);
  }
  else
  {
    die "Tier 1 or Annotation files were missing from $data_dir!\n";
  }

  $stats{'tier1_snvs_not_included'} = 0 if(!$stats{'tier1_snvs_not_included'});

  $stats{'tier1_snvs'} = 0 if(!$stats{'tier1_snvs'});
  $stats{'tier1_snvs_not_included'} = 0 if(!$stats{'tier1_snvs_not_included'});
  $stats{'tier1_snvs_with_annotation'} = 0 if(!$stats{'tier1_snvs_with_annotation'});
  $stats{'tier1_snvs_written'} = 0 if(!$stats{'tier1_snvs_written'});
  $stats{'tier1_indels'} = 0 if(!$stats{'tier1_indels'});
  $stats{'tier1_indels_not_included'} = 0 if(!$stats{'tier1_indels_not_included'});
  $stats{'tier1_indels_with_annotation'} = 0 if(!$stats{'tier1_indels_with_annotation'});
  $stats{'tier1_indels_written'} = 0 if(!$stats{'tier1_indels_written'});

  print $stats{'tier1_snvs'} . " tier 1 SNVs\n";
  print $stats{'tier1_snvs_not_included'} . " not included in target list\n";
  print $stats{'tier1_snvs_with_annotation'} . " met criteria and had annotation\n";
  print $stats{'tier1_snvs_written'} . " were written to MAF file\n\n";

  print $stats{'tier1_indels'} . " tier 1 Indels\n";
  print $stats{'tier1_indels_not_included'} . " not included in target list\n";
  print $stats{'tier1_indels_with_annotation'} . " met criteria and had annotation\n";
  print $stats{'tier1_indels_written'} . " were written to MAF file\n";

}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_mutations
{
  my $variant_file = shift(@_);
  my $input = new FileHandle ($variant_file);
  my $lineCounter = 0;

  my %mutations = ();

  while (<$input>)
  {
    chomp;
    my $line = $_;
    $lineCounter++;

    (my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);

    $var = iupac_to_base($ref, $var);

    my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);

    $mutations{$key} = $line;
  }

  close($input);


  return(%mutations);
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_annotation
{
  my $variant_file = shift(@_);
  my $input = new FileHandle ($variant_file);
  my $lineCounter = 0;

  my %annotation = ();

  while (<$input>)
  {
    chomp;
    my $line = $_;
    $lineCounter++;

    (my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);

    my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);

    $annotation{$key} = $line;
  }

  close($input);


  return(%annotation);
}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_dbsnps
{
  my $variant_file = shift(@_);
  my $input = new FileHandle ($variant_file);
  my $lineCounter = 0;

  print "Parsing $variant_file\n";

  my %mutations = ();

  while (<$input>)
  {
    chomp;
    my $line = $_;
    $lineCounter++;

    (my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
    my @lineContents = split(/\t/, $line);
    my $numContents = @lineContents;

    my $dbsnp_rs_id = $lineContents[$numContents - 1];
    my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);

    $mutations{$key} = $dbsnp_rs_id;
  }

  close($input);

  print "$lineCounter dbSNPs loaded\n";

  return(%mutations);
}



#############################################################
# IUPAC to base - convert IUPAC code to variant base
#
#############################################################

sub iupac_to_base
{
  (my $allele1, my $allele2) = @_;

  return($allele2) if($allele2 eq "A" || $allele2 eq "C" || $allele2 eq "G" || $allele2 eq "T");

  if($allele2 eq "M")
  {
    return("C") if($allele1 eq "A");
    return("A") if($allele1 eq "C");
  }
  elsif($allele2 eq "R")
  {
    return("G") if($allele1 eq "A");
    return("A") if($allele1 eq "G");
  }
  elsif($allele2 eq "W")
  {
    return("T") if($allele1 eq "A");
    return("A") if($allele1 eq "T");
  }
  elsif($allele2 eq "S")
  {
    return("C") if($allele1 eq "G");
    return("G") if($allele1 eq "C");
  }
  elsif($allele2 eq "Y")
  {
    return("C") if($allele1 eq "T");
    return("T") if($allele1 eq "C");
  }
  elsif($allele2 eq "K")
  {
    return("G") if($allele1 eq "T");
    return("T") if($allele1 eq "G");
  }

  return($allele2);
}


#############################################################
# ParseBlocks - takes input file and parses it
#
#############################################################

sub iupac_to_genotype
{
  (my $allele1, my $allele2) = @_;

  return($allele2 . $allele2) if($allele2 eq "A" || $allele2 eq "C" || $allele2 eq "G" || $allele2 eq "T");

  if($allele2 eq "M")
  {
    return("AC") if($allele1 eq "A");
    return("CA") if($allele1 eq "C");
  }
  elsif($allele2 eq "R")
  {
    return("AG") if($allele1 eq "A");
    return("GA") if($allele1 eq "G");
  }
  elsif($allele2 eq "W")
  {
    return("AT") if($allele1 eq "A");
    return("TA") if($allele1 eq "T");
  }
  elsif($allele2 eq "S")
  {
    return("GC") if($allele1 eq "G");
    return("CG") if($allele1 eq "C");
  }
  elsif($allele2 eq "Y")
  {
    return("TC") if($allele1 eq "T");
    return("CT") if($allele1 eq "C");
  }
  elsif($allele2 eq "K")
  {
    return("TG") if($allele1 eq "T");
    return("GT") if($allele1 eq "G");
  }

  return($allele2);
}

#############################################################
# ParseBlocks - takes input file and parses it
#
#############################################################
sub trv_to_mutation_type
{
  my $trv_type = shift(@_);

  return("Missense_Mutation") if($trv_type eq "missense");
  return("Nonsense_Mutation") if($trv_type eq "nonsense");
  return("Nonstop_Mutation") if($trv_type eq "nonstop");
  return("Silent") if($trv_type eq "silent");
  return("Splice_Site") if($trv_type eq "splice_site" || $trv_type eq "splice_site_del" || $trv_type eq "splice_site_ins");
  return("Frame_Shift_Del") if($trv_type eq "frame_shift_del");
  return("Frame_Shift_Ins") if($trv_type eq "frame_shift_ins");
  return("In_Frame_Del") if($trv_type eq "in_frame_del");
  return("In_Frame_Ins") if($trv_type eq "in_frame_ins");
  return("RNA") if($trv_type eq "rna");
  return("3'UTR") if($trv_type eq "3_prime_untranslated_region");
  return("5'UTR") if($trv_type eq "5_prime_untranslated_region");
  return("3'Flank") if($trv_type eq "3_prime_flanking_region");
  return("5'Flank") if($trv_type eq "5_prime_flanking_region");

  warn "Unknown mutation type $trv_type\n";
  return("Unknown");
}

################################################################################################
# Execute - the main program logic
#
################################################################################################

sub byChrPos
{
  my ($chrom_a, $pos_a) = split(/\t/, $a);
  my ($chrom_b, $pos_b) = split(/\t/, $b);

  $chrom_a =~ s/X/23/;
  $chrom_a =~ s/Y/24/;
  $chrom_a =~ s/MT/25/;
  $chrom_a =~ s/[^0-9]//g;

  $chrom_b =~ s/X/23/;
  $chrom_b =~ s/Y/24/;
  $chrom_b =~ s/MT/25/;
  $chrom_b =~ s/[^0-9]//g;

  $chrom_a <=> $chrom_a
  or
  $pos_a <=> $pos_b;
}

1;
