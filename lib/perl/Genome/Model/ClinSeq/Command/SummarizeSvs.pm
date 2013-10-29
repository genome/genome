package Genome::Model::ClinSeq::Command::SummarizeSvs;

#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Genome::Model::ClinSeq;
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::ClinSeq::Command::SummarizeSvs {
    is => 'Command::V2',
    has_input => [
        builds => { 
            is => 'Genome::Model::Build::SomaticVariation',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'somatic variation build(s) to summarize SVs from',
        },
        cancer_annotation_db => {
            is => 'Genome::Db::Tgi::CancerAnnotation',
            example_values => [$Genome::Model::ClinSeq::DEFAULT_CANCER_ANNOTATION_DB_ID],
            doc => 'cancer annotation data',
        },
        outdir => { 
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written', 
        },
    ],
    has_output => [
        fusion_output_file => {
            is => 'Text',
            doc => 'Path to output file',
        },
    ],
    doc => 'summarize the SVs of somatic variation build',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq summarize-svs --outdir=/tmp/  119390903 --cancer-annotation-db

genome model clin-seq summarize-svs --outdir=/tmp/  id=119390903

genome model clin-seq summarize-svs --outdir=/tmp/  model.id=2882504846

genome model clin-seq summarize-svs --outdir=/tmp/  "model.name='ALL1/AML103 - tumor/normal - wgs somatic variation - build37/hg19 - nov 2011 PP'"

genome model clin-seq summarize-svs --outdir=/tmp/  'id in [119390903,119517041]'

EOS
}

sub help_detail {
    return <<EOS
Summarize structural variants for one or more somatic variation builds 

(put more content here)
EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (-e $self->outdir && -d $self->outdir) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['outdir'],
	                                          desc => "Outdir: " . $self->outdir . " not found or not a directory",
                                          );
  }
  return @errors;
}

sub execute {
  my $self = shift;
  my @builds = $self->builds;
  my $outdir = $self->outdir;

  unless ($outdir =~ /\/$/){
    $outdir .= "/";
  }

  #Get Entrez and Ensembl data for gene name mappings
  #TODO: Create official versions of these data on allocated disk
  #Directory of gene lists for various purposes
  my $cancer_annotation_db = $self->cancer_annotation_db;
  my $clinseq_annotations_dir = $cancer_annotation_db->data_directory; 
  my $gene_symbol_lists_dir = $clinseq_annotations_dir . "/GeneSymbolLists/";
  my $entrez_ensembl_data = &loadEntrezEnsemblData(-cancer_db => $cancer_annotation_db);
  my $symbol_list_names = &importSymbolListNames('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-verbose'=>0);
  my $master_list = $symbol_list_names->{master_list};
  my @symbol_list_names = sort {$master_list->{$a}->{order} <=> $master_list->{$b}->{order}} keys %{$master_list};
  my $gene_symbol_lists = &importGeneSymbolLists('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-symbol_list_names'=>\@symbol_list_names, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);

  my $somatic_build_count = scalar(@builds);
  for my $somatic_build (@builds) {

    #If there is more than one somatic variation build supplied... create sub-directories for each
    my $build_outdir;
    if ($somatic_build_count > 1){
      $build_outdir = $outdir . $somatic_build->id . "/";
      mkdir ($build_outdir);
    }else{
      $build_outdir = $outdir;
    }

    my $build_dir = $somatic_build->data_directory;
    my $tumor_bam = $somatic_build->tumor_bam;
    my $normal_bam = $somatic_build->normal_bam;

    #Initially all of the following will be based on BreakDancear/SquareDancer output in the SomaticVariation output
    #e.g. /gscmnt/gc8002/info/model_data/2882679738/build119517041/variants/sv/union-union-sv_breakdancer_1.2__5-sv_breakdancer_1.2__6-sv_squaredancer_0.1__4/
    #squaredancer.svs.merge.file.annot

    #Make a copy of the annotated somatic SVs file and place in the outdir
    my $fusion_candidate_outfile = $build_outdir . "CandidateSvCodingFusions.tsv";
    my $sv_annot_search1 = $build_dir . "/variants/sv/union-union*/svs.hq.merge.annot.somatic";
    my $sv_annot_search2 = $build_dir . "/variants/sv/union-sv*/svs.hq.merge.annot.somatic";
     
    my $sv_annot_file = 'NULL';
    
    my $sv_annot_file1 = `ls $sv_annot_search1 2>/dev/null` || "NULL";
    chomp($sv_annot_file1);

    my $sv_annot_file2 = `ls $sv_annot_search2 2>/dev/null` || "NULL";
    chomp($sv_annot_file2);

    if (-e $sv_annot_file1){
      $sv_annot_file = $sv_annot_file1;
    }elsif (-e $sv_annot_file2){
      $sv_annot_file = $sv_annot_file2;
    }

    #Produce a simplified list of SVs gene fusion pairs (e.g. BCR-ABL1) - where type is fusion, and ORF affecting
    my %data;
    if (-e $sv_annot_file){

      #Copy the SV annot file to the clinseq working dir
      my $new_sv_annot_file = $build_outdir . "svs.hq.merge.annot.somatic";
      my $cp_cmd = "cp $sv_annot_file $build_outdir";
      unless (-e $new_sv_annot_file){
        Genome::Sys->shellcmd(cmd => $cp_cmd);
      }

      open (SV_ANNO, "$sv_annot_file") || die "\n\nCould not open SV annotation file: $sv_annot_file\n\n";
      my $l = 0;
      while(<SV_ANNO>){
        chomp($_);
        my @line = split("\t", $_);
        #Grab the somatic 'CTX' events, that are candidate 'Fusion' and 'AffectCoding'
        if ($_ =~ /\s+CTX\s+/ && $_ =~ /Fusion/ && $_ =~ /AffectCoding/){
          my $annot_string = $line[15];
          my $coord_string = $line[18];
          my $gene_pair = '';
          my $gene1 = '';
          my $gene2 = '';
          if ($annot_string =~ /^Gene\:(\S+?)\|(\S+?)\,/){
            $gene_pair = "$1-$2";
            $gene1 = $1;
            $gene2 = $2;
          }
          my $transcript1 = '';
          my $transcript2 = '';
          if ($annot_string =~ /\,(\w{2}\_\d+)\|.*\-(\w{2}\_\d+)\|/){
            $transcript1 = $1;
            $transcript2 = $2;
          }
          $l++;
          my @coords = split(",", $coord_string);
          my $mapped_gene_name1 = &fixGeneName('-gene'=>$gene1, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
          my $mapped_gene_name2 = &fixGeneName('-gene'=>$gene2, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
          my $record = "$gene_pair\t$gene1\t$gene2\t$coords[0]\t$coords[1]\t$mapped_gene_name1\t$mapped_gene_name2";
          $data{$l}{record} = $record;
          $data{$l}{coords1} = $coords[0];
          $data{$l}{coords2} = $coords[1];
          $data{$l}{gene1} = $gene1;
          $data{$l}{gene2} = $gene2;
          $data{$l}{mapped_gene1} = $mapped_gene_name1;
          $data{$l}{mapped_gene2} = $mapped_gene_name2;
          $data{$l}{transcript1} = $transcript1;
          $data{$l}{transcript2} = $transcript2;

        }
      }
    }else{
      $self->status_message("Could not find: SV annotation file by:\n$sv_annot_search1\n$sv_annot_search2");
    }

    #Annotate the genes of this file to help identify genes of interest (e.g. kinases, etc.)...
    foreach my $l (keys %data){
      my $gene1_name = $data{$l}{gene1};
      my $gene2_name = $data{$l}{gene2};

      foreach my $gene_symbol_type (keys %{$gene_symbol_lists}){
        my $gene_symbols = $gene_symbol_lists->{$gene_symbol_type}->{symbols};
        $data{$l}{$gene_symbol_type} = 0;
        if ($gene_symbols->{$gene1_name}){
          $data{$l}{$gene_symbol_type}++;
        }
        if ($gene_symbols->{$gene2_name}){
          $data{$l}{$gene_symbol_type}++;
        }
      }
    }

    #TODO: Use the coordinates of each fusion to produce pairoscope plots showing the support for each rearrangement
    #Relevant pairoscope option:
    #Usage:   pairoscope [options] <align.bam> <chr> <start> <end> <align2.bam> <chr2> <start2> <end2> 
    # -q INT    minimum mapping quality [0]
    # -p FLAG   output in pdf instead of png
    # -b INT    size of buffer around region to include [100]
    # -n FLAG   print the normal pairs too
    # -o STRING filename of the output png
    # -W INT    Width of the document [1024]
    # -H INT    Height of the document [768]
    # -g STRING bam file of exons for gene models
    # -u INT    upper bound of the insert size for a normal read [2147483647]
    # -l INT    lower bound of the insert size for a normal read[-1]
    # -m INT    minimum size of an event to display. Translocations are always displayed[0]
    # -P FLAG   Only display reads with both mates mapped in the graph
    # -t STRING list of transcripts to display

    #pairoscope -P -m 100000 -b 10000 /gscmnt/gc4074/info/model_data/2875512598/build111201309/alignments/111286100.bam 6 28884413 28884413 /gscmnt/gc4074/info/model_data/2875512598/build111201309/alignments/111286100.bam 17 28298036 28298036 -o HG1_TRIM27-EFCAB5_DNA.png
    
    #Set up the outdir and a temp file that will be used to gather read support counts dumped by pairoscope
    my $pairoscope_outdir = $build_outdir . "pairoscope/";
    mkdir($pairoscope_outdir);
    my $pairoscope_tmp_tumor_file = $pairoscope_outdir . "pairoscope_tumor.tmp";
    my $pairoscope_tmp_normal_file = $pairoscope_outdir . "pairoscope_normal.tmp";

    #Pairoscope run parameters
    my $min_mapping_quality = 1;
    my $min_event_size = 100000;
    my $flank = 10000;
    my $offset = 1000; #To make it easier to see the arcs...
    my $params_string = "-P -q $min_mapping_quality -m $min_event_size -b $flank";
    foreach my $l (keys %data){    
      my $gene1 = $data{$l}{gene1};
      my $gene2 = $data{$l}{gene2};
      my $coords1 = $data{$l}{coords1};
      my $coords2 = $data{$l}{coords2};
      my $transcript1 = $data{$l}{transcript1};
      my $transcript2 = $data{$l}{transcript2};

      #TODO: Figure out with Dave L. how this transcript option works.  Need to specify an Exons BAM with -g option?
      #TODO: Where do the annotations in the sv.annot file in somatic variation results come from
      my $transcript_string = '';
      if ($transcript1 && $transcript2){
        $transcript_string = "-t $transcript1,$transcript2";
      }

      my ($chr1, $start1, $end1, $chr2, $start2, $end2);
      if ($coords1 =~ /chr(\w+)\:(\d+)\-(\d+)/){
        $chr1 = $1;
        $start1 = ($2-$flank)-$offset;
        $end1 = ($3+$flank)-$offset;
      }

      if ($coords2 =~ /chr(\w+)\:(\d+)\-(\d+)/){
        $chr2 = $1;
        $start2 = ($2-$flank)+$offset;
        $end2 = ($3+$flank)+$offset;
      }

      #Pairoscope plot for Tumor BAM
      my $tumor_outfile = $pairoscope_outdir . "$gene1-$gene2"."_tumor_$l".".png";
      my $pairoscope_tumor_cmd = "pairoscope $params_string $transcript_string -o $tumor_outfile $tumor_bam $chr1 $start1 $end1 $tumor_bam $chr2 $start2 $end2 2>$pairoscope_tmp_tumor_file";
      Genome::Sys->shellcmd(cmd => $pairoscope_tumor_cmd);
      my $pairoscope_tumor_reads = 0;
      if (-e $pairoscope_tmp_tumor_file){
        open (TMP1, "$pairoscope_tmp_tumor_file") || die "\n\nCould not open tmp file: $pairoscope_tmp_tumor_file\n\n";
        while(<TMP1>){
          chomp($_);
          my @line = split("\t", $_);
          next unless (scalar(@line) == 4);
          $pairoscope_tumor_reads++;
        }
        close (TMP1);
      }
      $data{$l}{pairoscope_tumor_reads} = $pairoscope_tumor_reads;
      #print "\n\nTUMOR: $gene1 $coords2 $transcript1\t\t$gene2 $coords2 $transcript2\ttumor count = $pairoscope_tumor_reads\n$pairoscope_tumor_cmd";

      #Pairoscope plot for Normal BAM
      my $normal_outfile = $pairoscope_outdir . "$gene1-$gene2"."_normal_$l".".png";
      my $pairoscope_normal_cmd = "pairoscope $params_string $transcript_string -o $normal_outfile $normal_bam $chr1 $start1 $end1 $normal_bam $chr2 $start2 $end2 2>$pairoscope_tmp_normal_file";
      Genome::Sys->shellcmd(cmd => $pairoscope_normal_cmd);
      my $pairoscope_normal_reads = 0;
      if (-e $pairoscope_tmp_normal_file){
        open (TMP2, "$pairoscope_tmp_normal_file") || die "\n\nCould not open tmp file: $pairoscope_tmp_normal_file\n\n";
        while(<TMP2>){
          chomp($_);
          my @line = split("\t", $_);
          next unless (scalar(@line) == 4);
          $pairoscope_normal_reads++;
        }
        close (TMP2);
      }
      $data{$l}{pairoscope_normal_reads} = $pairoscope_normal_reads;
      #print "\n\nNORMAL: $gene1 $coords2 $transcript1\t\t$gene2 $coords2 $transcript2\tnormal count = $pairoscope_normal_reads\n$pairoscope_normal_cmd";

      #Clean-up the pairoscope temp files
      Genome::Sys->shellcmd(cmd => "rm -f $pairoscope_tmp_tumor_file");
      Genome::Sys->shellcmd(cmd => "rm -f $pairoscope_tmp_normal_file");
    }

    #Print out a new file containing the extra annotation columns
    open (FUSION_OUT, ">$fusion_candidate_outfile") || die "\n\nCould not open fusion outfile\n\n";
    my @gene_symbol_list_names = sort {$gene_symbol_lists->{$a}->{order} <=> $gene_symbol_lists->{$b}->{order}} keys %{$gene_symbol_lists};
    my $gene_symbol_list_name_string = join("\t", @gene_symbol_list_names);
    my $header_line = "gene_pair\tgene1\tgene2\tcoord1\tcoord2\tmapped_gene_name1\tmapped_gene_name2\tpairoscope_tumor_reads\tpairoscope_normal_reads";
    print FUSION_OUT "$header_line\t$gene_symbol_list_name_string\n";
    foreach my $l (sort {$a <=> $b} keys %data){
      my @tmp;
      foreach my $gene_symbol_list_name (@gene_symbol_list_names){
        push (@tmp, $data{$l}{$gene_symbol_list_name});
      }
      my $new_cols_string = join("\t", @tmp);
      print FUSION_OUT "$data{$l}{record}\t$data{$l}{pairoscope_tumor_reads}\t$data{$l}{pairoscope_normal_reads}\t$new_cols_string\n";
    }
    close(FUSION_OUT);
    $self->fusion_output_file($fusion_candidate_outfile);


    #TODO: Create a Stats.tsv file summarizing basic statistics of the sv annotations file
    
    #TODO: Identify the deletion regions.  Create a display for each deletion showing coverage across that region in tumor and normal

    #TODO: How reliable are our SVs.  Is there a way to identify a 'high confidence' set? 

    #TODO: Apply additional annotation strategies to the somatic SVs identified. 
    #- What are the genes/transcripts affected by the breakpoints of deletions, inversions, translocations?



  }
  $self->status_message("\n\n");

  return 1;
}

1;


