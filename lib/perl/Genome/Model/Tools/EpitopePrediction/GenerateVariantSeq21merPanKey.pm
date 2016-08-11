package Genome::Model::Tools::EpitopePrediction::GenerateVariantSeq21merPanKey;

use strict;
use warnings;
use Genome;

my $DEFAULT_TRV_TYPE = 'missense';


class Genome::Model::Tools::EpitopePrediction::GenerateVariantSeq21merPanKey {
    is => 'Genome::Model::Tools::EpitopePrediction::Base',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The input file is a tab-separated (TSV) output file from gmt annotate variant-protein. For more info, gmt annotate variant-protein --help',
        },
        output_file => {
            is => 'Text',
            doc => 'The output FASTA file to write 21mer sequences for wildtype(WT) and mutant(MT) proteins',
        },
        
        key_file => {
            is => 'Text',
            doc => 'The output lookup key for wildtype(WT) and mutant(MT) fasta headers',
        },
        trv_type => {
            is => 'Text',
            is_optional => 1,
            doc => 'The type of mutation you want to output eg missense,nonsense',
            # the current code only works on missense. Will need furthur development for other trv_types.
            default_value => $DEFAULT_TRV_TYPE,
        },
        
    ],
};


sub help_brief {
    "FOR NETHMHCPAN : Outputs a FASTA file for the wildtype(WT) and mutant(MT) proteins 21-mer sequences for MHC Class I epitope prediction. Assigns each protein to a Random Key and outputs the lookup in a separate key file. This is useful when the epitope prediction software truncates the FASTA header (eg. NetMHCPan)",
}



sub execute {
    my $self = shift;
    #my $tmp_dir = Genome::Sys->create_temp_directory();
    
    my $key_fh = Genome::Sys->open_file_for_writing($self->key_file);   
    my $input_fh = Genome::Sys->open_file_for_reading($self->input_file);
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    
 
 	my $i=1;
    while (my $line = $input_fh->getline) {
		chomp $line;
		$line =~ s/[*]$//g;
		my @prot_arr =  split(/\t/, $line);
		if ( $prot_arr[15] =~ /^p.([A-Z])(\d+)([A-Z])/ && $prot_arr[13] eq $self->trv_type)
		{
	#		open(OUT, '>>', $self->output_file) or die $!;
			my $wt_aa = $1;
			my $position = ($2 - 1);
			my $mt_aa = $3;
			my $wt_seq = $prot_arr[21];
			my @arr_wt_seq = split('',$wt_seq);

    		if ($1 ne $arr_wt_seq[$position])
    		{
    		next;
    		#TO DO :print OUT $prot_arr[0]."\t".$prot_arr[1]."\t".$prot_arr[2]."\t".$prot_arr[6]."\t".$1."\t".$2."\t".$3."\t".$prot_arr[11]."\t".$arr_wt_seq[$position]."\n";
    		}
   
   		 	if ($1 eq $arr_wt_seq[$position])
    		{
    		 my @mt_arr;
    	 	 my @wt_arr;
         		if ($position < 10)
        		{
           			@wt_arr = @arr_wt_seq[ 0 ... 20];
           			$arr_wt_seq[$position]=$mt_aa;
           			@mt_arr = @arr_wt_seq[ 0 ... 20];
		     		print $key_fh "WT_".$i."\t"."WT.".$prot_arr[6].".".$wt_aa.($position+1).$mt_aa."\n";
		     		print $output_fh ">WT_".$i."\n";
           			print $output_fh ( join "", @wt_arr);
           			print $output_fh "\n";
           			print $key_fh "MT_".$i."\t"."MT.".$prot_arr[6].".".$wt_aa.($position+1).$mt_aa."\n";
           			print $output_fh ">MT_".$i."\n";
           			print $output_fh ( join "", @mt_arr);
           			print $output_fh "\n";
        		}  
        		elsif ($position > ($#arr_wt_seq -10))
        		{	
           			@wt_arr = @arr_wt_seq[ $#arr_wt_seq -21 ... $#arr_wt_seq];
           			$arr_wt_seq[$position]=$mt_aa;
           			@mt_arr = @arr_wt_seq[ $#arr_wt_seq -21 ... $#arr_wt_seq];
           			print $key_fh "WT_".$i."\t"."WT.".$prot_arr[6].".".$wt_aa.($position+1).$mt_aa."\n";
		     		print $output_fh ">WT_".$i."\n";
           			print $output_fh ( join "", @wt_arr);
           			print $output_fh "\n";
           			print $key_fh "MT_".$i."\t"."MT.".$prot_arr[6].".".$wt_aa.($position+1).$mt_aa."\n";
           			print $output_fh ">MT_".$i."\n";
           			print $output_fh ( join "", @mt_arr);
           			print $output_fh "\n";
        		}	
       		 	elsif (($position >= 10) && ($position  <= ($#arr_wt_seq -10)))
        		{
           			@wt_arr = @arr_wt_seq[ $position-10 ... $position+10];
           			$arr_wt_seq[$position]=$mt_aa;
           			@mt_arr = @arr_wt_seq[ $position-10 ... $position+10];
           			print $key_fh "WT_".$i."\t"."WT.".$prot_arr[6].".".$wt_aa.($position+1).$mt_aa."\n";
		     		print $output_fh ">WT_".$i."\n";
           			print $output_fh ( join "", @wt_arr);
           			print $output_fh "\n";
           			print $key_fh "MT_".$i."\t". "MT.".$prot_arr[6].".".$wt_aa.($position+1).$mt_aa."\n";
           			print $output_fh ">MT_".$i."\n";
           			print $output_fh ( join "", @mt_arr);
           			print $output_fh "\n";
       				}
        		else 
        		{
        		print $output_fh "NULL"."\t".$position."\n";
        		}
        	#	OUT->close;
        		
    }
} 
$i++; 
    }
1;

}
return 1;



__END__
1	chromosome_name
2	start
3	stop
4	reference
5	variant
6	type
7	gene_name
8	transcript_name
9	transcript_species
10	transcript_source
11	transcript_version
12	strand
13	transcript_status
14	trv_type
15	c_position
16	amino_acid_change
17	ucsc_cons
18	domain
19	all_domains
20	deletion_substructures
21	transcript_error
22	wildtype_amino_acid_sequence

1 11
2 17801103
3 17801103	
4 G	
5 C	
6 SNP	
7 KCNC1	
8 NM_001112741.1	
9 human	
10 genbank	
11 58_37c	
12 +1	
13 reviewed	
14 missense	
15 c.1605	
16 p.P553L	
17 0.992	
18 NULL	
19 HMMSmart_SM00225,HMMPfam_K_tetra,HMMPfam_Ion_trans,superfamily_POZ domain,superfamily_Voltage-gated potassium channels	
20 -	
21 no_errors	
22 MGQGDESERIVINVGGTRHQTYRSTLRTLPGTRLAWLAEPDAHSHFDYDPRADEFFFDRHPGVFAHILNYYRTGKLHCPADVCGPLYEEELAFWGIDETDVEPCCWMTYRQHRDAEEALDSFGGAPLDNSADDADADGPGDSGDGEDELEMTKRLALSDSPDGRPGGFWRRWQPRIWALFEDPYSS








