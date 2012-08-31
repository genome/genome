package Genome::Model::Tools::Bed::Convert::Indel::StrelkaToBed;

use strict;
use warnings;

use Genome;
use Genome::Info::IUB; 

class Genome::Model::Tools::Bed::Convert::Indel::StrelkaToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Indel'],
    has_param => [
        limit_variants_to => { is => 'Text', valid_values => ['hq','lq'], is_optional => 1,
                              doc => 'set to "hq" or "lq" to only get variants which pass or fail filter, respectively' },
    ],
};

sub help_brief {
    "Tools to convert strelka Indel format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert indel strelka-to-bed --source snps_all_sequences --output snps_all_sequences.bed
EOS
}

sub help_detail {                           
    return <<EOS
    This is a small tool to take Indel calls in strelka format and convert them to a common BED format (using the first five columns).
EOS
}

sub process_source {
    my $self = shift;
    
    my $input_fh = $self->_input_fh;
    my $limit_variants_to = $self->limit_variants_to;
    
    #Skip comments and check header line and die if we find a problem
    while(my $line = <$input_fh>) {
      chomp($line);
      if ($line =~ /^\#\#/){
        next();
      }elsif($line =~ /^\#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+NORMAL\s+TUMOR/){
        last();
      }else{
        die "Bad header: $line";
      }
    }

    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR

    #Insertion example
    #1  54710 . A AT  . QSI_ref IC=6;IHP=8;NT=ref;QSI=20;QSI_NT=20;RC=5;RU=T;SGT=ref->het;SOMATIC;TQSI=1;TQSI_NT=1  DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50  22:22:20,21:0,0:2,2:18.64:0.00:0.00 56:56:38,38:12,15:6,9:51.5:0.00:0.00

    #Deletion example
    #1  404661  . CGT C . QSI_ref;Repeat  IC=22;IHP=2;NT=ref;QSI=10;QSI_NT=10;RC=23;RU=GT;SGT=ref->het;SOMATIC;TQSI=1;TQSI_NT=1 DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50  10:10:10,11:0,0:0,1:7.97:0.00:0.00  49:49:32,33:9,11:8,8:34.25:0.00:0.00

    #Continue processing the data lines
    while(my $line = <$input_fh>) {
      chomp($line);
      my @fields = split("\t", $line);
        my ($chromosome, $position, $id, $reference_base, $variant_base, $quality, $filter, $info, $format, $normal_info, $tumor_info) = @fields;
        if (index ($format, "DP:") != 0){
          die "Unrecognized format string on line: $line";
        }
        if ($limit_variants_to) {
            if ($limit_variants_to eq 'hq' and $filter ne 'PASS') {
                next;
            }
            elsif ($limit_variants_to eq 'lq' and $filter eq 'PASS') {
                next;
            }
        }
        my $somatic_score;
        if ($info =~ /\bQSI\=(\d+)\b/){
          $somatic_score = $1;
        }else{
          die "Can not parse somatic score from info field from line: $line";
        }
        my $tumor_depth;
        if ($tumor_info =~ /^(\d+)\:/){
          $tumor_depth = $1;
        }else{
          die "Can not parse tumor depth from tumor info field from line: $line";
        }

        #We want the tumor *genotype* using ambiguiety codes
        #This can probably be inferred from the $reference_base and $variant_base fields in the Strelka results above
        #e.g. A   AT  would be */T  (for insertions the start and end coordinates are always equal) 
        #e.g. CGT C   would be GT/* (for deletions the end coordinate is greater than the start coordinate by the size of the deletion)
        #It seems like either the reference or variant base is always 1 base.  
        #For insertions, the reference base is size 1.  For deletions, the variant base is size 1.
        #      
        
        my $ref_size = length($reference_base);
        my $var_size = length($variant_base);
        my $start_position = $position;

        my $inserted_bases;
        my $deleted_bases;
        my $end_position = 0;
        if ($ref_size == 1){
          #Insertion
          $deleted_bases = "*";
          $inserted_bases = substr($variant_base, 1);
          $end_position = $start_position;
        }elsif ($var_size == 1){
          #Deletion
          $deleted_bases = substr($reference_base, 1);
          $inserted_bases = "*";
          $end_position = $start_position + length($deleted_bases);
        }else{
          die "Expected either reference or variant string to be size 1 but it was not. Line: $line";
        }
        
        #print "\n\n$line\n$reference_base\t$variant_base\t$indel_string";


        #position => 1-based position of the Indel
        #BED uses 0-based position of and after the event
        $self->write_bed_line($chromosome, $position - 1, $end_position - 1, $deleted_bases, $inserted_bases, $somatic_score, $tumor_depth);

     }
    
    return 1;
}


1;
