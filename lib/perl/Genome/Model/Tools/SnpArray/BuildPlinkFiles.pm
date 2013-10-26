package Genome::Model::Tools::SnpArray::BuildPlinkFiles;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ModelGroup - Build Genome Models for Capture Datasets
#
#  AUTHOR:    Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#  CREATED:  12/09/2009 by D.K.
#  CHANGE LOG:
#  06/06/2011 [CK] Added support for new ref-align directory structure to access snp files (variants/snvs.hq.bed)
#
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;

use Genome; # using the namespace authorizes Class::Autouse to lazy-load modules under it
use Genome::Model::Tools::Analysis::Helpers qw(
    code_to_genotype
);

## Declare global statistics hash ##
my %stats = ();

class Genome::Model::Tools::SnpArray::BuildPlinkFiles {
  is => 'Command',

  has => [ # specify the command's single-value properties (parameters) <---
    group_id    => { is => 'Text', doc => "ID of model group" , is_optional => 0},
    map_file    => { is => 'Text', doc => "Location of map file in chrom, id, cM, position format" , is_optional => 0},
    sample_status_file    => { is => 'Text', doc => "Two column tab-delimited file of sample name and case(2) or control (1) status" , is_optional => 1},
    output_file    => { is => 'Text', doc => "Output file in PED format" , is_optional => 0},    
  ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Build PLINK files from SNP array data for a ref align model group"
}

sub help_synopsis {
    return <<EOS
Build PLINK files from SNP array data for a ref align model group
EXAMPLE:  gmt snp-array build-plink-files --group-id 661
EOS
}

sub help_detail { # this is what the user will see with the longer version of help. <---
    return <<EOS

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {
     my $self = shift;
     
     my $group_id = $self->group_id;
     my $map_file = $self->map_file;
     
     ## Keep stats in a single hash ##
     
     my %stats = ();
     
     ## Load sample statuses ##
     
     my %sample_status = ();
     print "Loading sample statuses...\n";
     %sample_status = load_sample_status_file($self->sample_status_file) if($self->sample_status_file);
     my $num_cases = my $num_controls = 0;
     foreach my $sample (keys %sample_status)
     {
          my $status = $sample_status{$sample};
          if($status eq "1")
          {
               $num_controls++;
          }
          elsif($status eq "2")
          {
               $num_cases++;
          }
     }
     
     print "$num_cases cases, $num_controls controls in sample set\n";
     
     ## Load the map file ##
     print "Loading map file...\n";
     my @map = load_map_file($map_file);
     
     ## Save the positions of all map keys ##
     
     my %map_keys = ();
     foreach my $snp (@map)
     {
          $stats{'num_map_snps'}++;
          my ($chrom, $id, $cm, $position) = split(/\t/, $snp);
          my $key = join("\t", $chrom, $position);
          $map_keys{$key} = $id;
     }
     
     print $stats{'num_map_snps'} . " SNPs in map file\n";

     open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
     open(NOTUSED, ">" . $self->output_file . ".notused") or die "Can't open outfile: $!\n";

     my %sample_snp_keys = ();
     
     ## Get the models in each model group ##
     
     my $model_group = Genome::ModelGroup->get($group_id);
     my @models = $model_group->models;
     
     foreach my $model (@models)
     {
          $stats{'models_in_group'}++;
          print $stats{'models_in_group'} .' '. $model->id . "\n";
          
          my $model_id = $model->genome_model_id;
          my $subject_name = $model->subject_name;
          my $genotype_model = $model->genotype_microarray_model;
          next if not $genotype_model;
          my $genotype_build = $genotype_model->last_succeeded_build;
          next if not $genotype_build;
          my $genotype_build_dir = $genotype_build->data_directory;
          
          ## PRoceed if we have a build directory ##
          my $genotype_file = $genotype_build_dir . "/gold_snp.v2.bed";
          $stats{'models_with_genotypes'}++;
          if(-e $genotype_file)  
          {
            
              print "Loading genotypes for $subject_name...\n";
              my %sample_genotypes = load_genotype_file($genotype_file);

              my $gts_present = my $gts_missing = 0;
              my $subject_status = $sample_status{$subject_name};
              $subject_status = 0 if(!$subject_status);
              my $sample_plink = join("\t", $subject_name, $subject_name, 0, 0, 1, $subject_status);
              $stats{'models_with_genotypes_status_' . $subject_status}++;
              #Go through map positions
              my %genotype_used = ();
              foreach my $snp (@map)
              {
                  my ($chrom, $id, $cm, $position) = split(/\t/, $snp);
                  my $key = join("\t", $chrom, $position);
                  if($sample_genotypes{$key})
                  {
                      $genotype_used{$key} = 1;
                      $gts_present++;
                      ## Separate out alleles ##

                      my $allele1 = substr($sample_genotypes{$key}, 0, 1);
                      my $allele2 = substr($sample_genotypes{$key}, 1, 1);
#                              my ($allele1, $allele2) = split(/\//, $sample_genotypes{$key});
                      $sample_plink .= "\t$allele1 $allele2";
                  }
                  else
                  {
                      $gts_missing++;
                      $sample_plink .= "\t0 0";
                  }

              }

              print OUTFILE $sample_plink . "\n";

              ## Go through sample SNPs and find map positions that are missing ##

              foreach my $key (keys %sample_genotypes)
              {
                  $sample_snp_keys{$key}++;
                  if(!$genotype_used{$key})
                  {
                      print NOTUSED "$key\n";
                  }
              }

              print join("\t", $subject_name, $genotype_file, $gts_present) . "\n";
#                    exit(0) if($stats{'models_with_genotypes'} >= 4);
          }
      }

      close(OUTFILE);

      print $stats{'models_in_group'} . " models in group\n" if($stats{'models_in_group'});
      print $stats{'models_with_genotypes'} . " had succeeded genotype builds\n";
      print $stats{'models_with_genotypes_status_2'} . " cases\n";
      print $stats{'models_with_genotypes_status_1'} . " controls\n";
      print $stats{'models_with_genotypes_status_0'} . " unknowns\n";
      my $num_snps = my $num_map = my $num_novel = 0;

#    Uncomment this if we need to build a map file 
      #
#     open(OUTFILE, ">Sample-SNP-Keys.tsv");
#     foreach my $key (keys %sample_snp_keys)
#     {
#          $num_snps++;
#          if($map_keys{$key})
#          {
#               $num_map++;
#               print OUTFILE join("\t", $key, "MAP") . "\n";
#          }
#          else
#          {
#               $num_novel++;
#               print OUTFILE join("\t", $key, "NOVEL") . "\n";
#          }
#     }
#     close(OUTFILE);

      print "$num_snps SNPs seen in samples\n";
      print "$num_map were in the map file\n";
      print "$num_novel were not\n";

      return 1;

  }





#############################################################
# load the map file
  #
#############################################################

  sub load_map_file
  {
      my $FileName = shift(@_);

      my @snps = ();
      my $snpCounter = 0;

      my $input = new FileHandle ($FileName);
      my $lineCounter = 0;


      while (<$input>)
      {
          chomp;
          my $line = $_;
          $lineCounter++;
          my ($chrom, $id, $cm, $position) = split(/\t/, $line);
          $snps[$snpCounter] = join("\t", $chrom, $id, $cm, $position);
          $snpCounter++;
      }
      close($input);

      return(@snps);
  }




#############################################################
# load the map file
  #
#############################################################

  sub load_sample_status_file
  {
      my $FileName = shift(@_);

      my %status = ();

      my $input = new FileHandle ($FileName);
      my $lineCounter = 0;

      while (<$input>)
      {
          chomp;
          my $line = $_;
          my ($sample_name, $sample_status) = split(/\t/, $line);
          $status{$sample_name} = $sample_status;
      }

      close($input);

      return(%status);
  }


#############################################################
# load the map file
  #
#############################################################

  sub load_genotype_file
  {
      my $FileName = shift(@_);

      my %genotypes = ();

      my $input = new FileHandle ($FileName);
      my $lineCounter = 0;

      while (<$input>)
      {
          chomp;
          my $line = $_;
          $lineCounter++;
          my ($chrom, $chr_start, $chr_stop, $alleles) = split(/\t/, $line);
          $chrom = 23 if ($chrom eq "X");
          $chrom = 24 if($chrom eq "Y");
          $chrom = 26 if($chrom eq "MT");
          my $key = join("\t", $chrom, $chr_stop);
          my ($ref, $cns) = split(/\//, $alleles);
          my $genotype = code_to_genotype($cns);

          $genotypes{$key} = $genotype;
      }

      close($input);

      return(%genotypes);
  }

1;
