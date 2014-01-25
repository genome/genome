#!/usr/bin/env genome-perl
use strict;
use warnings;
use Getopt::Long;
use GSCApp;
use FileHandle;
use above 'Genome';

use Bio::Seq;
use Bio::SeqIO;

if(@ARGV<1) {die "$0  annotation_file\n";}

#MPSampleData::DBI::myinit("dbi:Oracle:dwrac","mguser_prd"); #switch to production by default

#my $file=new FileHandle("/gscmnt/sata148/info/medseq/llin/LOH/BC1/sequencing/BRC1_unique_mutations.csv");
my $file=new FileHandle($ARGV[0]);
#my $file=new FileHandle("/gscmnt/sata423/info/medseq/llin/aml/copy_number_analysis/all_Jackie/Copy_of_AML_Gold1_2_Recurrence_manual_review_combined_wo_dbsnps_confirmantion_review_100430.050210.csv");
open OUT, ">$ARGV[0].sift.polyphen";
while(<$file>)
{
    chomp;
    my $line=$_;
    #if (/^chr/i)
    if (/chromosome/)
    {
	print OUT "$line\tpolyphen\tsift\n";
	#print "$line\n";
	next;
    }
    my @lines=split /\t/,$line;
    #my $chromosome=$lines[0];
    #my $start=$lines[1];
    #my $end=$lines[2];
    #my $type=$lines[5];
    #my $transcript=$lines[7];
    #my $aac=$lines[14];
    #my $polyphen=$line[24];
    #my $sift=$line[25];
    #my ($chromosome,$start,$end,$ref,$var,$type,$gene,$transcript,$tr_source,$tr_version,$strand,$tr_status,$trv_type,$c_pos,$aac,$ucsc_cons,$domain,$normal_freq,$tumor_freq,$met_freq,$xeno_freq,$primary_call,$met_call,$xeno_call,$polyphen,$sift,@others)=split/\t/;
    
    my ($chromosome,$start,$end,$ref,$var,$type,$gene,$transcript_name,$tr_species,$tr_source,$tr_version,$strand,$tr_status,$trv_type,$c_pos,$aac,$ucsc_cons,$domain,$all_domains,$deletion_substructures,$transcript_error)=split/\t/;

    #my ($chromosome,$start,$stop,$ref,$var,$type,$somatic_status,$comment,$man_ref_status,$gene,$transcript_name,$tr_species,$tr_scource,$tr_version,$strand,$tr_status,$trv_type,$c_position,$aac,$ucsc_cons,$domain,$all_domains,$deletion_substructures,$transcript_error)=split /\t/, $line;

    unless ($trv_type =~ /missense/i) {
        print OUT "$line\tNULL\tNULL\n";
        next;
    }
    #next if (length($aac)==0);
    #$chromosome=5;
    #$start=135422725;
    #$end=$start;
    #$aac="p.E576K";
    #$transcript="NM_000358";
print "transcript: $transcript_name\n";

    my $transcript_iterator = Genome::Model->get(name => "NCBI-human.combined-annotation")->build_by_version("58_37c_v2")->transcript_iterator(where => [chrom_name => $chromosome, start_position => $start]);

    my $transcript;
    my $count;
    until ($transcript){
        my $t = $transcript_iterator->next;
        unless ($t){
            die "transcript iterator exhausted after $count transcripts.  $transcript_name not found!";
        }
        $count++;
        $transcript = $t if $t->transcript_name eq $transcript_name;
    }

    my $prediction=get_prediction($transcript,$aac,1,1);
    my $polyphen=$prediction->{"polyphen"};
    my $sift=$prediction->{"sift"};

    #print "$chromosome\t$start\t$end\t$ref\t$var\t$type\t$gene\t$transcript_name\t$tr_species\t$tr_source\t$tr_version\t$strand\t$tr_status\t$trv_type\t$c_pos\t$aac\t$ucsc_cons\t$domain\t$all_domains\t$polyphen\t$sift\n";
    print OUT "$line\t$polyphen\t$sift\n";

}


sub get_prediction{
    my ($transcript,$pro_var,$polyphen,$sift)=@_;
    my $predict={"polyphen"=>"NULL","sift"=>"NULL"};
    return $predict unless($pro_var =~ /p\.[A-Z]\d+[A-Z]/ );
    return $predict if($pro_var =~ /[B|J|O|U|X|*|Z]/ );
    $pro_var =~ /p\.([A-Z])(\d+)([A-Z])/;
    my $a1=$1;
    my $num=$2;
    my $a2=$3;

    my $protein = $transcript->protein;
    unless ($protein){
        print "no protein found for transcript ".$transcript->transcript_name."\n";
    }
    my $pr_name=$protein->protein_name;
# write out the protein sequence file
    my $out = Bio::SeqIO->new(-file => ">$pr_name.tmp_proseq.fasta" , '-format' => 'fasta');
    my $seq = Bio::Seq->new(-seq=>$protein->amino_acid_seq,-display_id => $protein->protein_name);
    $out->write_seq($seq);

# prepare substitution file
    open(SUBSTITUTION,">$pr_name.tmp_pph.substitution") or die "Can't open $pr_name.tmp_pph.substitution. $!";
    print SUBSTITUTION "0\t0\t",$protein->protein_name,"\t$num\t$a1\t$a2\n";
    close(SUBSTITUTION);
    ` awk '{print \$5\$4\$6 }' $pr_name.tmp_pph.substitution > $pr_name.tmp_sift.substitution `;

# run polyphen
    if($polyphen){
        my $pphdire=qq(/gscuser/ndees/scripts/PolyPhen.1.13.sh);
        my $polycmd = "$pphdire -s $pr_name.tmp_proseq.fasta  $pr_name.tmp_pph.substitution";
        #print "$polycmd \nRunning...\n\n";
        `$polycmd > $pr_name.tmp_pphprediction.out`;
# update polyphen result
        #print "uploading polyphen results....\n";
        open(POLY, "<$pr_name.tmp_pphprediction.out") or goto label;
        while (<POLY>) {last;}
        while (<POLY>) {
            chomp;
            if($_=~/Error:/) {$predict->{"polyphen"}="Error"; next;}
            my @fields = split("\t");
            $fields[6]=~tr/ /_/;
            $predict->{"polyphen"}= $fields[6];
        }
        close(POLY);
#print "finish loading polyphen results....\n";
    }

    if($sift){
#  my $db_location = "/gscmnt/sata180/info/medseq/xshi/nr"; 
#    my $db_location = "/gsc/var/lib/blastdb/ncbi/nr";
        my $db_location = "/gscmnt/sata180/info/medseq/xshi/uniprot/hs_swall";
#  my $db_location = "/gscmnt/sata180/info/medseq/xshi/uniprot/swissprot";
        my $siftcmd = "SIFT.csh  $pr_name.tmp_proseq.fasta $db_location $pr_name.tmp_sift.substitution";
#print "\n$siftcmd\nRunning...\n\n";
        `$siftcmd > $pr_name.tmp_siftprediction.out`;
# update sift result to database
#print "uploading sift results....\n";
        open(SIFT, "<$pr_name.tmp_proseq.SIFTprediction") or goto label;
# while (<SIFT>) {last;}
        while (<SIFT>) {
            chomp;
            my @fields = split("\t");
            next unless($fields[0]=~/[A-Z]\d+[A-Z]/);
            $fields[2]="" unless(defined $fields[2]);
            $predict->{"sift"}=lc($fields[1])."<>".$fields[2];
        }

        close(SIFT);
#print "finish loading sift results....\n";
    }


    label: #print "finishing prediction\n"; 
    `rm -rf $pr_name.tmp_pphprediction.out $pr_name.tmp_siftprediction.out $pr_name.tmp_proseq.* $pr_name.tmp_pph.substitution $pr_name.tmp_sift.substitution`;
    return $predict;
}
