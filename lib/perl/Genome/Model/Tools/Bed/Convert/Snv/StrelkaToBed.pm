package Genome::Model::Tools::Bed::Convert::Snv::StrelkaToBed;

use strict;
use warnings;

use Genome;
use Genome::Info::IUB; 

class Genome::Model::Tools::Bed::Convert::Snv::StrelkaToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Snv'],
    has_param => [
        limit_variants_to => { is => 'Text', valid_values => ['hq','lq'], is_optional => 1,
                              doc => 'set to "hq" or "lq" to only get variants which pass or fail filter, respectively' },
    ],
};

sub help_brief {
    "Tools to convert strelka SNV format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert snv strelka-to-bed --source snps_all_sequences --output snps_all_sequences.bed
EOS
}

sub help_detail {                           
    return <<EOS
    This is a small tool to take SNV calls in strelka format and convert them to a common BED format (using the first five columns).
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
    #1       10231   .       C       A       .       QSS_ref NT=ref;QSS=1;QSS_NT=1;SGT=AC->AC;SOMATIC;TQSS=2;TQSS_NT=2       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    32:4:8:0:0,3:28,60:0,0:0,1     84:6:69:0:7,21:71,192:0,0:0,1

    #Continue processing the data lines
    while(my $line = <$input_fh>) {
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
        if ($info =~ /\bQSS\=(\d+)\b/){
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
        #This can probably be obtained from the 'SGT' field in the Strelka results above
        #e.g. AA->AG where 'AG' is the predicted heterozygous somatic genotype and would need to be converted to the correct IUB code
        my @info_fields = split(";", $info);
        my $sgt = $info_fields[3];
        my $genotype_string;
        my $ambig_code;
        if ($sgt =~ /SGT\=(\w+)\-\>(\w+)/){
          my $nor_genotype = $1;
          my $tum_genotype = $2;
          $ambig_code = Genome::Info::IUB->string_to_iub($tum_genotype);
          $ambig_code || die "Failed to obtain ambiguity code for tumor genotype: $tum_genotype from line: $line";
        }else{
          die "Unable to process the 4th field of an info block, expected SGT field from line: $line";
        }

        #position => 1-based position of the SNV
        #BED uses 0-based position of and after the event
        $self->write_bed_line($chromosome, $position - 1, $position, $reference_base, $ambig_code, $somatic_score, $tumor_depth);
     }
    
    return 1;
}


1;
