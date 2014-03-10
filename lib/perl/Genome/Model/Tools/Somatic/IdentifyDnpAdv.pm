package Genome::Model::Tools::Somatic::IdentifyDnpAdv;

use strict;
use warnings;

use Genome;
use Command;
use Genome::Info::IUB;
use IO::File;
use POSIX;
use Sort::Naturally;
use List::Util qw(max);

my $DEFAULT_VERSION = '0.2';
my $READCOUNT_COMMAND = 'bam-readcount';

class Genome::Model::Tools::Somatic::IdentifyDnpAdv {
    is => 'Command',
    has => [
    bam_readcount_version => {
        is => 'Version',
        is_input=>1, 
        is_optional => 1,
        default_value => $DEFAULT_VERSION,
        doc => "Version of bam-readcount to use, default is $DEFAULT_VERSION"
    },
    bam_readcount_params => {
        is => 'String',
        is_optional => 1,
        is_input=>1, 
        default_value => "-q 1",
        doc => "Parameters to pass to bam-readcount",
    },
    snp_input_file => {
        type => 'String',
        is_optional => 0,
        is_input=>1,
        doc => 'List of sites in input file 1_based file format to look for DNPs. This must be sorted by chromosome and coordinate.',
        default => '',
    },
    proportion => {
        type => 'Float',
        is_optional => 1,
        default => 0.1,
        doc => 'Proportion of reads supporting the DNP required for the site to be considered a DNP',
    },
    bam_file => {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'File from which to retrieve reads. Must be indexed.',
    },
    'min_mapping_quality' => {
            type => 'String',
            default => '40',
            is_optional => 1,
            is_input => 1,
            doc => 'minimum average mapping quality threshold for high confidence call',
    },
    'min_somatic_score' => {
            type => 'String',
            default => '40',
            is_optional => 1,
            is_input => 1,
            doc => 'minimum somatic quality threshold for high confidence call',
     },
     anno_lc_file =>{
        type => 'String',
        is_optional => 0,
        is_output=>1,
        doc => '1_based format low confidence sites. This must be sorted by chromosome and coordinate.',
        default => '',
    },
    bed_lc_file => {
        type => 'String',
        is_optional => 0,
        is_output=>1,
        doc => '0_based format low confidence sites. This must be sorted by chromosome and coordinate.',
        default => '',
    },
    anno_hc_file => {
        type => 'String',
        is_optional => 0,
        is_output=>1,
        doc => '1_based high confidence sites. This must be sorted by chromosome and coordinate.',
        default => '',
    },
    bed_hc_file => {
        type => 'String',
        is_optional => 0,
        is_output=>1,
        doc => '0_based high confidence sites. This must be sorted by chromosome and coordinate.',
        default => '',
    },
    ]
};

my %READCOUNT_VERSIONS = (
    '0.2' => $ENV{GENOME_SW} . '/samtools/readcount/readcount-v0.2/' . $READCOUNT_COMMAND,
);

sub help_brief {
    "Scans an annotation file, finds adjacent sites and then identifies if these are DNPs."
}

sub help_detail {
    <<'HELP';
This is a simple script which operated by identifying adjacent sites in a SORTED annotation file (really, just chr, start, stop, ref, cns, and type are need in that order). And then it breadk the IUB code and check all the possible recombination and determine if the alleles are linked in the same reads. If so, wrapping them into a single DNP and TNP event. The other sites are simply printed as is. Currently, only DNPs and TNPs are examined. Any continus SNPs sites > 3 will still be handled as SNPs. The input file was intended to be the pre-annotation "adapted file" of SNVs from the somatic pipeline file and, therefore, should NOT contain indels.The output is bed_format, that can be put into fast-tiering  
HELP
}

#Code should operate as follows
#Scan the snp file
#Upon finding adjacent snps
#retrieve reads for the two snps
#Check and see if they are in cis ie in the same read together

sub execute {
    my $self=shift;
    my $snp_file = $self->snp_input_file;

    #check on architecture
    unless (POSIX::uname =~ /64/) {
        $self->error_message("This script requires a 64-bit system to run samtools");
        die;
    }

    #check on BAM file
    unless(-e $self->bam_file && !-z $self->bam_file) {
        $self->error_message($self->bam_file . " does not exist or is of zero size");
        die;
    }

    unless(-e $self->bam_file . ".bai") {
        $self->error_message("Tumor bam must be indexed");
        die;
    }
    
    #check on input snp file
    my $fh = IO::File->new($snp_file, "r");
    unless($fh) {
        $self->error_message("Couldn't open $snp_file: $!"); 
        die;
    }

    my ($last_chr,$last_stop,$last_pos, $last_ref, $last_cns, $last_type,$last_score); 
    my ($chr,$stop, $pos, $ref, $cns, $type,$score);
    my @rest=();
    my %SNP;  
    
    #the following logic assumes that you have a single position per line
    # looking for candidate NNP sites 
    while(my $line = $fh->getline) {
      chomp $line;
      ($chr, $pos,$stop, $ref, $cns, $type, $score, @rest) = split /\t/, $line;
      if ($last_chr && $last_pos){
          if ($last_chr eq $chr){
              if (($pos-$last_stop)==1){
                  $last_stop=$stop;
                  $last_ref=$last_ref.$ref;
                  $last_cns=$last_cns.$cns;
                  $last_score=$last_score.",".$score;
                  $last_type++;  
              }else{              
                  $SNP{$last_chr}{$last_pos}{$last_stop}{$last_ref}{$last_cns}{$last_type}=$last_score;
                  $last_chr = $chr;
                  $last_pos = $pos;
                  $last_stop= $stop;
                  $last_ref = $ref;
                  $last_cns = $cns;
                  $last_score=$score;
                  $last_type=1;
              } 
          }else{
              $SNP{$last_chr}{$last_pos}{$last_stop}{$last_ref}{$last_cns}{$last_type}=$last_score;
              $last_chr = $chr;
              $last_pos = $pos;
              $last_stop= $stop;
              $last_ref = $ref;
              $last_cns = $cns;
              $last_score=$score;
              $last_type=1;
           }
       }else{
              $last_chr = $chr;
              $last_pos = $pos;
              $last_stop= $stop;
              $last_ref = $ref;
              $last_cns = $cns;
              $last_score=$score;
              $last_type=1;
       } 
   }
   $SNP{$last_chr}{$last_pos}{$last_stop}{$last_ref}{$last_cns}{$last_type}=$score;
   $fh->close;

   #sort candidate variants according to Chromosome and then Start
    my $n=0; my @candidate=();
    for my $CHR (nsort keys %SNP){
        for my $POS (sort {$a <=>$b} keys %{$SNP{$CHR}}){
            for my $STOP (sort { $a cmp $b } keys %{$SNP{$CHR}{$POS}}){
                   for  my $REF (sort { $a cmp $b } keys %{$SNP{$CHR}{$POS}{$STOP}}){
                       for  my $CNS (sort { $a cmp $b } keys %{$SNP{$CHR}{$POS}{$STOP}{$REF}}){
                              for my $TYPE (sort {$a <=> $b} keys %{$SNP{$CHR}{$POS}{$STOP}{$REF}{$CNS}}){
                                  my $SCORE= $SNP{$CHR}{$POS}{$STOP}{$REF}{$CNS}{$TYPE};
                                  if ($TYPE >1){
                                      push @candidate,"$CHR\t$POS\t$STOP\t$REF\t$CNS\t$TYPE\t$SCORE";
                                      $n++;
                                  }else{
                                      push @candidate, "$CHR\t$POS\t$STOP\t$REF\t$CNS\t$TYPE\t$SCORE";
                                  }
                              }
                        }
                   }
              }
         }
    }
    $self->debug_message("nNNP candidate number: $n");
    my ($t_begin, $t_end);
    $t_begin=time();
    $self->debug_message("Begin find_dnp_from_candidate: $t_begin");
  
    my @result=$self->find_dnp_from_candidate(\@candidate);

    $t_end=time();
    $self->debug_message("Finish find_dnp_from_candidate: $t_end");

    # open FILEHANDLE for output bed format and anno format
    my $bedfile_hc = IO::File->new($self->bed_hc_file, "w");
    unless($bedfile_hc) {
        $self->error_message("Unable to open " . $self->bed_hc_file . " for writing. $!");
        die;
    }
    my $annofile_hc = IO::File->new($self->anno_hc_file, "w");
    unless($annofile_hc) {
        $self->error_message("Unable to open " . $self->anno_hc_file . " for writing. $!");
        die;
    }
    my $bedfile_lc = IO::File->new($self->bed_lc_file, "w");
    unless($bedfile_lc) {
        $self->error_message("Unable to open " . $self->bed_lc_file . " for writing. $!");
        die;
    }
    my $annofile_lc = IO::File->new($self->anno_lc_file, "w");
    unless($annofile_lc) {
        $self->error_message("Unable to open " . $self->anno_lc_file . " for writing. $!");
        die;
    }   

    # Do bam-readcount to check the mapping quality of each bases
    $t_end=time();
    $self->debug_message("Begin make tmp readcount file: $t_end");
     
    my ($tfh,$temp_path) = Genome::Sys->create_temp_file;
    unless($tfh) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }
    $temp_path =~ s/\:/\\\:/g;

   
    #read sniper and skip indels
    my %all_lines;
    for my $line (@result) {
        chomp $line;
        my ($chr, $start, $stop, $ref, $var, $type, $somatic_score, @annotation_columns) = split /\t/, $line;
        if (! exists $all_lines{$chr}{$start}{$stop} ){
            $all_lines{$chr}{$start}{$stop}=$line;
        }
#        if ($self->prepend_chr) {
#            $chr = "chr$chr";
#            $chr =~ s/MT$/M/;
#        };
        next if $ref eq "*";
        print $tfh "$line\n";
    }
    $tfh->close;
    $t_end=time();
    $self->debug_message("Begin run readcount: $t_end");
    my %count_line;
    my $readcount_command=sprintf("%s %s -l %s %s|",$self->readcount_path, $self->bam_readcount_params, $temp_path,$self->bam_file);
    $self->debug_message("Running: $readcount_command");
    my $readcounts = IO::File->new("$readcount_command") or die "can't open the file $readcount_command due to $!";
    my @readcounts= $readcounts->getlines;
    $t_end=time();
    $self->debug_message("End read readcount: $t_end");

    my $total_readcount=scalar(@readcounts);
    for my $count_line (@readcounts){
        chomp $count_line;
        my ($chr, $pos, $ref, $depth, @base_stats) = split /\t/, $count_line;
        $count_line{$chr}{$pos}=$count_line;
     }
    unless($readcounts->close()) {
        $self->error_message("Error running $ENV{GENOME_SW}/samtools/readcount/readcount-v0.2/bam-readcount");
        die;
    }
    $self->debug_message("The number of total sites for readcount: $total_readcount");
    $t_end=time();
    $self->debug_message("Begin parsing readcount and print outfile: $t_end");
    my $total_result= scalar(@result);
    $self->debug_message("The number of total sites before readcount: $total_result");

    # separate high and low confidence sites
    # For BWA, if (somatic_score >=40 && mapping_quality>=40) is fullfilled in any base of DNP/TNP, the DNP/TNP is high confidence

    EACH:
    for my $result (@result){
        chomp $result;
        my ($chr, $start, $stop, $ref, $var, $type, $somatic_score, @annotation_columns) = split /\t/, $result;
        READCOUNT:
        my $base_index=0;
        my $current_start =$start+$base_index;
        while($current_start <= $stop){
            if (exists $count_line{$chr}{$current_start}){
                my ($c_chr, $c_pos, $c_ref, $c_depth, @base_stats) = split /\t/, $count_line{$chr}{$current_start};
                my %bases;
                for my $base_stat (@base_stats){
                    my ($base, $reads, $avg_mq, $avg_bq)=split/\:/,$base_stat;
                    next if ($base eq "=" || $base eq "N");
                    $bases{$base}=$avg_mq;
                }
                if ($type eq "SNP"){
                    if (exists $bases{$var}){
                        if ($somatic_score >= $self->min_somatic_score && $bases{$var} >= $self->min_mapping_quality ){
                             print $annofile_hc  "$chr\t$start\t$stop\t$ref\t$var\t$type\t$somatic_score\t$bases{$var}\n";
                             my $bedfile_start=$start-1;
                             print $bedfile_hc  "$chr\t$bedfile_start\t$stop\t$ref\t$var\t$type\t$somatic_score\t$bases{$var}\n";
                         }else{
                             print $annofile_lc  "$chr\t$start\t$stop\t$ref\t$var\t$type\t$somatic_score\t$bases{$var}\n";
                             my $bedfile_start=$start-1;
                             print $bedfile_lc  "$chr\t$bedfile_start\t$stop\t$ref\t$var\t$type\t$somatic_score\t$bases{$var}\n";

                        }
                        next EACH;
                    }else{
                        $self->debug_message("cannot find $result in READCOUNT");
                    }
                }else{
                    my @somatic_scores=split/\,/, $somatic_score;
                    my $current_var=substr($var,$base_index,1);
                    my $current_somatic_score=$somatic_scores[$base_index];
                    if (exists $bases{$current_var}){
                         if ($current_somatic_score >= $self->min_somatic_score && $bases{$current_var} >= $self->min_mapping_quality ){
                             print $annofile_hc  "$chr\t$start\t$stop\t$ref\t$var\t$type\t$current_somatic_score\t$bases{$current_var}\n";
                             my $bedfile_start=$start-1;
                             print $bedfile_hc  "$chr\t$bedfile_start\t$stop\t$ref\t$var\t$type\t$current_somatic_score\t$bases{$current_var}\n";
                             next EACH;
                         }else{
                             if ($current_start==($stop-1)){
                                 print $annofile_lc  "$chr\t$start\t$stop\t$ref\t$var\t$type\t$current_somatic_score\t$bases{$current_var}\n";
                                 my $bedfile_start=$start-1;
                                 print $bedfile_lc  "$chr\t$bedfile_start\t$stop\t$ref\t$var\t$type\t$current_somatic_score\t$bases{$current_var}\n";
                                 next EACH;
                             }
                        }
                    }else{
                        $self->debug_message("cannot find $result in READCOUNT");
                    }
                }
            }
            $base_index++;
            $current_start =$start+$base_index;

        }
    }
    $bedfile_hc->close; $annofile_hc->close; $bedfile_lc->close; $annofile_lc->close;
    $t_end=time();
    $self->debug_message("End parsing readcount and print outfile: $t_end");
    if (($t_end-$t_begin) > 0) {
        my $avg = $total_result/($t_end-$t_begin);
        $self->debug_message("Average speed is $avg sites per second");
    }
    return 1;
}


#In the find_dnp_from_candidate, the somatic score for SNP remains the same
#keep the somatic scores for each sites in DNP and TNP for later to separate high/low confidence sites
sub find_dnp_from_candidate{
    my $self=shift;
    my $ref_array=shift;
    my (@candidate)=@{$ref_array};
    my @result=();
    my $n=0;
    for my $candidate (@candidate){
        $candidate=~s/\n//;
        my (@dnp, @dnp1,@dnp2,@snp, @tnp) = ((),(),(),(),()); 
        my ($chr, $pos, $stop, $ref, $cns, $type, $score) = split /\t/, $candidate;
        if ($type==1) { 
            my @var1=Genome::Info::IUB::variant_alleles_for_iub($ref,$cns);
            for my $var1(@var1){
                if ($var1 ne $ref) {
                 push @result, "$chr\t$pos\t$stop\t$ref\t$var1\tSNP\t$score\n";
                }
            }
        }elsif ($type==2){
            $n++;
#            print "$candidate\n";
            my $ref1=substr($ref,0,1); 
            my $ref2=substr($ref,1,1);
            my $cns1=substr($cns,0,1);
            my $cns2=substr($cns,1,1);
            my @var1=Genome::Info::IUB::variant_alleles_for_iub($ref1,$cns1); 
            my @var2=Genome::Info::IUB::variant_alleles_for_iub($ref2,$cns2); 
            my $stop1=$pos+1;
            my ($score1, $score2)=split/\,/,$score;
            for my $var1 (@var1){
                for my $var2 (@var2){
                   if ($self->is_dnp($chr,$pos,$var1,$stop1,$var2)){
                       if ($ref ne $var1.$var2){
                           my $score_dnp=$score1+$score2;
                           my $dnp= "$chr\t$pos\t$stop1\t$ref\t$var1$var2\tDNP\t$score\n";
                           push @dnp, $dnp; 
                       }
                   }
               }
            } 
            if ( (scalar@dnp) > 0){
                for my $tmp_dnp (@dnp){
                    push @result, $tmp_dnp;
                }
            }
        }elsif ($type ==3 ){
            $n++; @dnp1=(); @dnp2=();
            my ($score1,$score2,$score3)=split/\,/,$score;
            my $ref1=substr($ref,0,1); 
            my $ref2=substr($ref,1,1); 
            my $ref3=substr($ref,2,1);
            my $cns1=substr($cns,0,1);
            my $cns2=substr($cns,1,1);
            my $cns3=substr($cns,2,1); 
            my @var1=Genome::Info::IUB::variant_alleles_for_iub($ref,$cns1); 
            my @var2=Genome::Info::IUB::variant_alleles_for_iub($ref,$cns2); 
            my @var3=Genome::Info::IUB::variant_alleles_for_iub($ref,$cns3);
            my $stop1=$pos+1; my $stop2=$pos+2; 
            for my $var1 (@var1){
                for my $var2 (@var2){
                   if ($self->is_dnp($chr,$pos,$var1,$stop1,$var2)){
                       if ( $ref1.$ref2 ne $var1.$var2){
                           my $dnp1= "$chr\t$pos\t$stop1\t$ref1$ref2\t$var1$var2\tDNP\t$score1,$score2\n";
                           push @dnp1,$dnp1;
                       }
                   }
               }
            } 
            for my $var2 (@var2){
                for my $var3 (@var3){
                   if ($self->is_dnp($chr,$stop1,$var2,$stop2,$var3)){
                       if ($ref2.$ref3 ne $var2.$var3){ 
                           my $dnp2="$chr\t$stop1\t$stop2\t$ref2$ref3\t$var2$var3\tDNP\t$score2,$score3\n"; 
                           push @dnp2,$dnp2;
                       }
                   }
               }
            }
            if((scalar@dnp1)>0 && (scalar@dnp2)>0){
                   #FIXME Assume double dnp is TNP if the middle nucleotide is same. 
                   #Need to remove the false positive case, such as AT and TG are on different reads
                   #TODO sort the @tnp,@dnp,@snp together before print out. It now print out SNP>DNP>TNP
                   for my $dnp1(@dnp1){
                       for my $dnp2(@dnp2){
                           my ($chr1,$pos1,$stop1,$ref1,$cns1,$type1)=$dnp1;
                           my ($chr2,$pos2,$stop2,$ref2, $cns2,$type2)=$dnp2;
                           my $ref1_1=substr($ref1,0,1);
                           my $cns1_1=substr($cns1,0,1);
                           my $cns1_2=substr($cns1,1,1);
                           my $cns2_1=substr($cns2,0,1);
                           my $cns2_2=substr($cns2,1,1);
                           if ($cns1_2 eq $cns2_1 && $ref1_1.$ref2 ne $cns1.$cns2_2){
                               my $tnp_var=$cns1_1.$cns1_2.$cns2_1;
                               my $tnp_ref=substr($ref1,0,1).substr($ref2,0,1).substr($ref2,1,1);
                               my $tnp="$chr\t$pos1\t$stop2\t$tnp_ref\t$tnp_var\tTNP\t$score\n";
                               push @tnp,$tnp; 
                           }else{
                               push @dnp,$dnp1;
                               push @dnp,$dnp2;
                           }
                       }
                   }
            }elsif ((scalar@dnp1)>0 || (scalar@dnp2)>0){ 
                # only one DNP + SNP survive
                if ((scalar@dnp1)>0){
                    @dnp=@dnp1;
                    my $snp="$chr\t$stop2\t$stop2\t$ref3\t$cns3\tSNP\t$score3\n";
                    push @snp, $snp;
                }elsif((scalar@dnp2)>0){ 
                    @dnp=@dnp2; 
                    my $snp="$chr\t$pos\t$pos\t$ref1\t$cns1\tSNP\t$score1\n";
                    push @snp, $snp;
                }
            }else{
                # only 3 SNPs independently exists
                my $snp="$chr\t$pos\t$pos\t$ref1\t$cns1\tSNP\t$score1\n";
                push @snp, $snp;
                $snp="$chr\t$stop1\t$stop1\t$ref2\t$cns2\tSNP\t$score2\n";
                push @snp, $snp;
                $snp="$chr\t$stop2\t$stop2\t$ref3\t$cns3\tSNP\t$score3\n"; 
                push @snp, $snp;
            }
            my $new=undef;
            if ((scalar @snp) >0 ){
                for $new (@snp){ push @result,$new; }
            }
            if ((scalar @dnp)>0 ){
                for $new (@dnp){ push @result, $new; }
            }
            if ((scalar @tnp)>0 ){
                for $new (@tnp){ push @result, $new; }
            }
        }elsif ($type > 3){ 
            # FIXME CANNOT handle SNPs >3, make it back to SNPs
            $self->error_message("Unable to process more than TNP\n$candidate\, make it back to SNPsn");
            my ($CHR,$START,$STOP,$REF,$VAR,$TYPE)=split/\t/,$candidate;
            my $current=$START; my $n=0;
            while ( ($current+$n) < $STOP){
                my $start =$current+$n;
                my $cns=substr($VAR,$n,1);
                my $ref=substr($VAR,$n,1);
                my @scores=split/\,/,$score;
                my @var1=Genome::Info::IUB::variant_alleles_for_iub($ref,$cns);
                my $index=0;
                for my $var1(@var1){
                    my $score1=$scores[$index];
                    if ($var1 ne $ref) {
                        push @result, "$chr\t$pos\t$stop\t$ref\t$var1\tSNP\t$score1\n";
                    }
                }
                $n++;
            }
        }
    }
    return @result;
}

#I can't say why this is just a wrapper for the other program
sub is_dnp {
    my ($self, $chr, $pos1, $base1, $pos2, $base2) = @_;
    return $self->_determine_dnp_from_bam_reads($self->bam_file,$chr,$pos1,$base1,$pos2,$base2);
}

#This grabs the reads overlapping the positions
#and checks to see if they contain both potential DNP bases
sub _determine_dnp_from_bam_reads {
    my ($self, $alignment_file, $chr, $pos1, $base1, $pos2, $base2) = @_;
    unless(open(SAMTOOLS, "samtools view $alignment_file $chr:$pos1-$pos2 |")) {
        $self->error_message("Unable to open pipe to samtools view");
        return;
    }
    my ($reads, $reads_supporting_dnp) = (0,0);
    while( <SAMTOOLS> ) {
        chomp;
        my ($qname, $flag, $rname, $pos_read, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual, $RG, $MF, @rest_of_fields) = split /\t/;
        next if($mapq == 0); #only count q1 and above

        my $offset1 = $self->_calculate_offset($pos1, $pos_read, $cigar);
        next unless defined $offset1; #skip deletions
        my $offset2 = $self->_calculate_offset($pos2, $pos_read, $cigar);
        next unless defined $offset2; #skip deletions
        $reads++;
        
        if(uc(substr($seq,$offset1,1)) eq uc($base1) && uc(substr($seq,$offset2,1)) eq uc($base2)) {
            $reads_supporting_dnp++;
        } 
    }
    unless(close(SAMTOOLS)) {
        $self->error_message("Error running samtools");
        return -1;
    }
    if($reads_supporting_dnp/$reads > $self->proportion) {
         return 1;
    }else{
        return 0;
    }
}

#this calculates the offset of a position into a seqeunce string based on the CIGAR string specifying the alignment
sub _calculate_offset { 
    my $self = shift;
    my $pos = shift;
    my $read_pos = shift;
    my $cigar = shift;
    my $current_offset=0;
    my $current_pos=$read_pos;
    my @ops = $cigar =~ m/([0-9]+)([MIDNSHP])/g; 
    OP:
    while(my ($cigar_len, $cigar_op) =  splice @ops, 0, 2 ) {
        my $new_offset;
        my $last_pos=$current_pos;
        if($cigar_op eq 'M') {
            $current_pos+=$cigar_len;
            $current_offset+=$cigar_len;
        }
        elsif($cigar_op eq 'I') {
            $current_offset+=$cigar_len;
        }
        elsif($cigar_op eq 'D') {
            $current_pos+=$cigar_len;

        }
        elsif($cigar_op eq 'N') {
            #this is the same as a deletion for returning a base from the read
            $current_pos += $cigar_len;
        }
        elsif($cigar_op eq 'S') {
            #soft clipping means the bases are in the read, but the position (I think) of the read starts at the first unclipped base
            #Functionally this is like an insertion at the beginning of the read
            $current_offset+=$cigar_len;
        }
        elsif($cigar_op eq 'H') {
            #hard clipping means the bases are not in the read and the position of the read starts at the first unclipped base
            #Shouldn't do anything in this case, but ignore it
        }
        else {
            die("CIGAR operation $cigar_op currently unsupported by this module");
        }
        if($pos < $current_pos && $pos >= $last_pos) {
            if($cigar_op eq 'M') {
                my $final_adjustment = $current_pos - $pos;
                return $current_offset - $final_adjustment;
            }
            else {
                return;
            }
        }
    }
    #position didn't cross the read
    return; 
}

sub readcount_path {
    my $self = $_[0];
    return $self->path_for_readcount_version($self->bam_readcount_version);
}

sub available_readcount_versions {
    my $self = shift;
    return keys %READCOUNT_VERSIONS;
}

sub path_for_readcount_version {
    my $class = shift;
    my $version = shift;

    if (defined $READCOUNT_VERSIONS{$version}) {
        return $READCOUNT_VERSIONS{$version};
    }
    die('No path for bam-readcount version '. $version);
}

sub default_readcount_version {
    die "default bam-readcount version: $DEFAULT_VERSION is not valid" unless $READCOUNT_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

