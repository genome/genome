#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#For a group of ClinSeq models, converge various types of Cufflinks expression values together
#For example, allow for merging at the level of isoforms and genes
#Make it as generic as possible
#Deal with SampleType.  i.e. get both 'tumor' and 'normal' files if available

#Input:
#A list of Clinseq builds, models, or a Clinseq model-group

#Parameters:
#1.) the name of the file to be joined. e.g. 'isoforms.merged.fpkm.expsort.tsv'
#2.) the join column name (unique ID).  e.g. 'tracking_id'
#3.) the data column name containing the data to be used for the matrix.  e.g. 'FPKM'
#4.) a list of annotation column names.  Values from these will be taken from only file and appended to the end of the resulting matrix file

#Sanity checks:
#Make sure each file found has the neccessary columns specified by the user
#Make sure each file found has the same number of IDs, only unique IDs, and all the same IDs as every other file
#Make sure all values are defined.  If there are empty cells allow an option for these to be converted to NAs

#Output:
#A single expression matrix file
#In the output expression file, name the expression columns according to: CommonName_SampleType_BuildID
#Sort by primary ID

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::OriginalScripts::Converge qw(:all);

my $build_ids = '';
my $model_ids = '';
my $model_group_id = '';
my $target_file_name = '';
my $expression_subdir = '';
my $join_column_name = '';
my $data_column_name = '';
my $annotation_column_names = '';
my $outfile = '';
my $verbose = 0;

GetOptions ('build_ids=s'=>\$build_ids, 'model_ids=s'=>\$model_ids, 'model_group_id=s'=>\$model_group_id, 
            'target_file_name=s'=>\$target_file_name, 'expression_subdir=s'=>\$expression_subdir,
            'join_column_name=s'=>\$join_column_name, 'data_column_name=s'=>\$data_column_name, 'annotation_column_names=s'=>\$annotation_column_names,
            'outfile=s'=>\$outfile, 'verbose=i'=>\$verbose);

my $usage=<<INFO;
  Example usage: 

  Gene-level using isoforms merged to each gene
  convergeCufflinksExpression.pl  --model_group_id='25134'  --target_file_name='isoforms.merged.fpkm.expsort.tsv'  --expression_subdir='isoforms_merged'  --join_column_name='tracking_id'  --data_column_name='FPKM'  --annotation_column_names='ensg_name,mapped_gene_name,locus'  --outfile=Cufflinks_GeneLevel_Malat1Mutants.tsv  --verbose=1

  Transcript-level using isoforms individually
  convergeCufflinksExpression.pl  --model_group_id='25134'  --target_file_name='isoforms.fpkm.expsort.tsv'  --expression_subdir='isoforms'  --join_column_name='tracking_id'  --data_column_name='FPKM'  --annotation_column_names='gene_id,mapped_gene_name,locus'  --outfile=Cufflinks_IsoformLevel_Malat1Mutants.tsv  --verbose=1

  Specify *one* of the following as input (each model/build should be a ClinSeq model)
  --build_ids                Comma separated list of specific build IDs
  --model_ids                Comma separated list of specific model IDs
  --model_group_id           A singe genome model group ID

  Combines Cufflinks expression results from a group of Clinseq models into a single report:
  --target_file_name         The files to be joined across multiple ClinSeq models
  --expression_subdir        The expression subdir of Clinseq to use: ('genes', 'isoforms', 'isoforms_merged')

  --join_column_name         The primary ID to be used for joining values (IDs must be unique and occur in all files to be joined)
  --data_column_name         The data column to be used to create an expression matrix across the samples of the ClinSeq models
  --annotation_column_names  Optional list of annotation columns to append to the end of each line in the matrix (value will be taken from the first file parsed)
                             The values in these columns are assumed to be constant across the files being merged!
  --outfile                  Path of the output file to be written
  --verbose                  More descriptive stdout messages

  Test Clinseq model groups:
  25307                      BRAF inhibitor resistant cell lines
  25134                      MALAT1 mutant vs. wild type BRCs
  30176                      LUC RNA-seq vs. RNA-cap
  31430                      ALS, spinal cord vs. occiptal cortex

INFO

unless (($build_ids || $model_ids || $model_group_id) && $target_file_name && $expression_subdir && $join_column_name && $data_column_name && $outfile){
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
  exit(1);
}

#Get the models/builds
if ($verbose){print BLUE, "\n\nGet genome models/builds for supplied list", RESET;}
my $models_builds;
if ($build_ids){
  my @build_ids = split(",", $build_ids);
  $models_builds = &getModelsBuilds('-builds'=>\@build_ids, '-verbose'=>$verbose);
}elsif($model_ids){
  my @model_ids = split(",", $model_ids);
  $models_builds = &getModelsBuilds('-models'=>\@model_ids, '-verbose'=>$verbose);
}elsif($model_group_id){
  $models_builds = &getModelsBuilds('-model_group_id'=>$model_group_id, '-verbose'=>$verbose);
}else{
  print RED, "\n\nCould not obtains models/builds - check input to convergeCufflinksExpression.pl\n\n", RESET;
  exit();
}

#Get the annotation columns specified by the user:
my @annotation_columns;
if ($annotation_column_names){
  @annotation_columns = split(",", $annotation_column_names);
}
my $annotation_columns_string = join("\t", @annotation_columns);

#Get the desired files from each model
my %files = %{&getCufflinksFiles('-models_builds'=>$models_builds, '-target_file_name'=>$target_file_name, '-expression_subdir'=>$expression_subdir)};

#Get the list of distinct data columns that will be written
my %data_list;
foreach my $fc (keys %files){
  my $column_name = $files{$fc}{column_name};
  $data_list{$column_name}=1;
}
my @data_list_sort = sort keys %data_list; 
my $data_list_string = join("\t", @data_list_sort);

#Parse each file and build a hash keyed on the join id and output column name.  Store target data as a value.  Store annotation values
if ($verbose){print BLUE, "\n\nParse all files found for data '$data_column_name' joining on '$join_column_name", RESET;}
my $exp = &parseCufflinksFiles('-files'=>\%files, '-data_column_name'=>$data_column_name, '-join_column_name'=>$join_column_name, '-annotation_columns'=>\@annotation_columns);

#Print the output file
if ($verbose){print BLUE, "\n\nPrinting output to: $outfile", RESET;}
open (OUT, ">$outfile") || die "\n\nCould not open output file for writing: $outfile\n\n";
my $header;
if ($annotation_column_names){
  $header = "$join_column_name\t$data_list_string\t$annotation_columns_string";
}else{
  $header = "$join_column_name\t$data_list_string";
}
print OUT "$header\n";
foreach my $id (sort keys %{$exp}){
  my @data;
  foreach my $data_col (@data_list_sort){
    push(@data, $exp->{$id}->{$data_col}->{data});
  }
  my $data_string = join("\t", @data);
  if ($annotation_column_names){
    my $annotation_string = $exp->{$id}->{annotation};
    print OUT "$id\t$data_string\t$annotation_string\n";
  }else{
    print OUT "$id\t$data_string\n";
  }
}
close(OUT);

if ($verbose){
  print "\n\n";
}else{
  print "\n\nPrinted output to: $outfile\n\n";
}

exit();


#######################################################################################################################
#Get the desired files from each model                                                                                #
#######################################################################################################################
sub getCufflinksFiles{
  my %args = @_;
  my $models_builds = $args{'-models_builds'};
  my $target_file_name = $args{'-target_file_name'};
  my $expression_subdir = $args{'-expression_subdir'};

  my %files;
  my $fc = 0;
  if ($verbose){print BLUE, "\n\nGet all Cufflinks files within these builds that match $target_file_name", RESET;}
  my %mb = %{$models_builds->{cases}};
  foreach my $c (keys %mb){
    my $b = $mb{$c}{build};
    my $m = $mb{$c}{model};
    my $data_directory = $b->data_directory;
    my $subject_name = $b->subject->name;
    my $subject_common_name = $b->subject->common_name;
    my $build_id = $b->id;

    unless ($m->normal_rnaseq_model) {
        die "model " . $m->__display_name__ . " does not have a normal rnaseq model!";
    }

    my ($normal_rnaseq_subject, $tumor_rnaseq_subject); 
    $normal_rnaseq_subject = $m->normal_rnaseq_model->subject->name || "NULL";
    $tumor_rnaseq_subject = $m->tumor_rnaseq_model->subject->name || "NULL";

    #If the subject name is not defined, die
    unless ($subject_name){
      print RED, "\n\nCould not determine subject name for build: $build_id\n\n", RESET;
      exit(1);
    }

    my $final_name = "Unknown";
    if ($subject_name){$final_name = $subject_name;}
    if ($subject_common_name){$final_name = $subject_common_name;}

    if ($verbose){print BLUE, "\n\t$final_name\t$build_id\t$data_directory", RESET;}

    #/gscmnt/gc8002/info/model_data/2881869913/build120828540/BRC18/rnaseq/tumor/cufflinks_absolute/isoforms_merged
    my $ls_cmd = "ls $data_directory/*/rnaseq/*/cufflinks_expression_absolute/$expression_subdir/* 2>/dev/null";
    my @result = `$ls_cmd`;
    chomp(@result);
    my @files;
    foreach my $result (@result){
      if ($result =~ /$target_file_name$/){
        push(@files, $result);
      }
    }

    #Get the common name and subtype (e.g. tumor/normal) from the path
    foreach my $file (@files){
      if ($file =~ /$data_directory\/(.*)\/rnaseq\/(.*)\/cufflinks.*/){
        $fc++;
        $files{$fc}{path} = $file;
        $files{$fc}{subtype} = $2;
        $files{$fc}{final_name} = $final_name;
        $files{$fc}{subject_name} = $subject_name;
        $files{$fc}{build_id} = $build_id;
        if ($file =~ /normal/){
          $files{$fc}{subject_name2} = $normal_rnaseq_subject;
        }
        if ($file =~ /tumor/){
          $files{$fc}{subject_name2} = $tumor_rnaseq_subject;
        }

      }else{
        print RED, "\n\nCould not obtain subtype from file path:\n$file\n\n", RESET;
        exit(1);
      }
    }
  }

  #Determine a suitable name for data columns in the output file.  It must be distinct accross the builds.  Try these combinations in order:
  # final_name
  # final_name + subtype
  # final_name + subtype + subject_name
  # final_name + subtype + subject_name + build_id
  my $file_count =  keys %files;
  unless ($file_count > 1){
    print RED, "\n\nFound $file_count files to join (need at least 2 ...)\n\n", RESET;
    exit(1);
  }
  my %finalnames;
  my %finalnames_subtypes;
  my %finalnames_subtypes_subjectnames;
  my %finalnames_subtypes_subjectnames_subjectnames2;
  my %finalnames_subtypes_subjectnames_subjectnames2_buildids;
  foreach my $fc (keys %files){
    my $final_name = $files{$fc}{final_name};
    my $subtype = $files{$fc}{subtype};
    my $subject_name = $files{$fc}{subject_name};
    my $subject_name2 = $files{$fc}{subject_name2};
    my $build_id = $files{$fc}{build_id};
    my $id1 = $final_name;
    my $id2 = $id1 . "_"."$subtype";
    my $id3 = $id2 . "_"."$subject_name";
    my $id4 = $id3 . "_"."$subject_name2";
    my $id5 = $id4 . "_"."$build_id";
    $finalnames{$id1}=1;
    $finalnames_subtypes{$id2}=1;
    $finalnames_subtypes_subjectnames{$id3}=1;
    $finalnames_subtypes_subjectnames_subjectnames2{$id4}=1;
    $finalnames_subtypes_subjectnames_subjectnames2_buildids{$id5}=1;
    $files{$fc}{id1} = $id1;
    $files{$fc}{id2} = $id2;
    $files{$fc}{id3} = $id3;
    $files{$fc}{id4} = $id4;
    $files{$fc}{id5} = $id5;
  }

  my $id_choice;
  if ($file_count == keys %finalnames){
    $id_choice = "id1";
  }elsif($file_count == keys %finalnames_subtypes){
    $id_choice = "id2";
  }elsif($file_count == keys %finalnames_subtypes_subjectnames){
    $id_choice = "id3";
  }elsif($file_count == keys %finalnames_subtypes_subjectnames_subjectnames2){
    $id_choice = "id4";
  }elsif($file_count == keys %finalnames_subtypes_subjectnames_subjectnames2_buildids){
    $id_choice = "id5";
  }else{
    print RED, "\n\nCould not generate a suitable set of distinct labels for the data files to be joined...\n\n", RESET;
    exit();
  }

  foreach my $fc (keys %files){
    $files{$fc}{column_name} = $files{$fc}{"$id_choice"};
  }

  return(\%files);
}


#######################################################################################################################
#parseCufflinksFiles                                                                                                  #
#######################################################################################################################
sub parseCufflinksFiles{
  my %args = @_;
  my %files = %{$args{'-files'}};
  my $data_column_name = $args{'-data_column_name'};
  my $join_column_name = $args{'-join_column_name'};
  my @annotation_columns = @{$args{'-annotation_columns'}};

  my %exp;
  my $first_file = 1;
  foreach my $fc (keys %files){
    my $path = $files{$fc}{path};
    my $column_name = $files{$fc}{column_name};
    if ($verbose){print BLUE, "\n\tProcess: $path ($column_name)", RESET;}
    my %columns;
    my $header = 1;
    my $line_count = 0;
    open (IN, "$path") || die "\n\nCould not open input file: $path\n\n";
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header){
        my $p = 0;
        foreach my $name (@line){
          $columns{$name}{p} = $p;
          $p++;
        }
        #Check for requested columns
        unless(defined ($columns{$join_column_name})){print RED, "\n\nCould not find required join column: $join_column_name\n\n", RESET; exit(1);}
        unless(defined ($columns{$data_column_name})){print RED, "\n\nCould not find required data column: $data_column_name\n\n", RESET; exit(1);}
        foreach my $annotation_column (@annotation_columns){unless(defined ($columns{$annotation_column})){print RED, "\n\nCould not find required annotation column: $annotation_column\n\n", RESET; exit(1);}}
        $header = 0;
        next();
      }
      $line_count++;
      my $id = $line[$columns{$join_column_name}{p}];
      my $data = $line[$columns{$data_column_name}{p}];

      #Unless this is the first file parsed, the id must already be defined but the id-column_name value must not be
      unless ($first_file){
        unless (defined($exp{$id})){
          print RED, "\n\nFound an ID in a file that was not present in the first file parsed. ID: $id\n\n", RESET;
          #exit(1);
        }
        if (defined($exp{$id}{$column_name})){
          print RED, "\n\nFound a duplicate ID: $id\n\n", RESET;
          exit(1);
        }
      }

      $exp{$id}{$column_name}{data} = $data;

      if ($annotation_column_names){
        my @annotation_values;
        foreach my $annotation_column (@annotation_columns){
          my $annotation = $line[$columns{$annotation_column}{p}];
          push(@annotation_values, $annotation);
        }
        my $annotation_string = join("\t", @annotation_values);
        $exp{$id}{annotation} = $annotation_string;
      }
    }
    close(IN);

    #Check that the total number of records stored equals the total number of data lines in the input file (i.e. there were no duplicate IDs)
    my $data_count = keys %exp;
    unless ($data_count == $line_count){
      print RED, "\n\nData count ($data_count) did not match the expected line count ($line_count)!\n\n", RESET;
      #exit(1);
    }
    $first_file = 0;
  }

  return(\%exp);
}





