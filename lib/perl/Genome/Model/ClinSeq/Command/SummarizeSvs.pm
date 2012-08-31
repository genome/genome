package Genome::Model::ClinSeq::Command::SummarizeSvs;

#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Term::ANSIColor qw(:constants);

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
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
    ],
    doc => 'summarize the SVs of somatic variation build',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq summarize-builds --outdir=/tmp/  119390903

genome model clin-seq summarize-builds --outdir=/tmp/  id=119390903

genome model clin-seq summarize-builds --outdir=/tmp/  model.id=2882504846

genome model clin-seq summarize-builds --outdir=/tmp/  "model.name='ALL1/AML103 - tumor/normal - wgs somatic variation - build37/hg19 - nov 2011 PP'"

genome model clin-seq summarize-builds --outdir=/tmp/  'id in [119390903,119517041]'

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
	                                          desc => RED . "Outdir: " . $self->outdir . " not found or not a directory" . RESET,
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

    #Initially all of the following will be based on BreakDancear/SquareDancer output in the SomaticVariation output
    #e.g. /gscmnt/gc8002/info/model_data/2882679738/build119517041/variants/sv/union-union-sv_breakdancer_1.2__5-sv_breakdancer_1.2__6-sv_squaredancer_0.1__4/
    #squaredancer.svs.merge.file.annot

    #Make a copy of the annotated somatic SVs file and place in the outdir
    my $fusion_candidate_outfile = $outdir . "CandidateSvCodingFusions.tsv";
    my $sv_annot_search = $build_dir . "/variants/sv/*/svs.hq.merge.annot.somatic";
    my $sv_annot_file = `ls $sv_annot_search` || "NULL";
    chomp($sv_annot_file);

    #Produce a simplified list of SVs gene fusion pairs (e.g. BCR-ABL1) - where type is fusion, and ORF affecting
    if (-e $sv_annot_file){

      #grep -w CTX /gscmnt/gc8002/info/model_data/2882504846/build119390903/variants/sv/union-union-sv_breakdancer_1.2__5-sv_breakdancer_1.2__6-sv_squaredancer_0.1__4/squaredancer.svs.merge.file.annot | grep Fusion | grep AffectCoding
      open (SV_ANNO, "$sv_annot_file") || die "\n\nCould not open SV annotation file: $sv_annot_file\n\n";
      open (FUSION_OUT, ">$fusion_candidate_outfile") || die "\n\nCould not open fusion outfile\n\n";
      print FUSION_OUT "gene_pair\tgene1\tgene2\tcoord1\tcoord2\n";
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
          my @coords = split(",", $coord_string);
          print FUSION_OUT "$gene_pair\t$gene1\t$gene2\t$coords[0]\t$coords[1]\n";
        }
      }
      close(FUSION_OUT);
    }else{
      $self->status_message("Could not find: SV annotation file: $sv_annot_search");
    }

    #Annotate the genes of this file to help identify genes of interest (e.g. kinases, etc.)...

    #Use the coordinates of each fusion to produce pairoscope plots showing the support for each rearrangement

    #Identify the deletion regions.  Create a display for each deletion showing coverage across that region in tumor and normal

    #How reliable are our SVs.  Is there a way to identify a 'high confidence' set? 

    #Apply additional annotation strategies to the somatic SVs identified. 
    #What are the genes/transcripts affected by the breakpoints of deletions, inversions, translocations?



  }
  $self->status_message("\n\n");

  return 1;
}

1;


