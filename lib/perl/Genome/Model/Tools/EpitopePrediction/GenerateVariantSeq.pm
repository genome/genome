package Genome::Model::Tools::EpitopePrediction::GenerateVariantSeq;

use strict;
use warnings;

use Genome;
use Workflow;

my $DEFAULT_TRV_TYPE = 'missense';


class Genome::Model::Tools::EpitopePrediction::GenerateVariantSeq {
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
    "FOR NETMHC : Outputs a FASTA file for  for wildtype(WT) and mutant(MT) proteins 21-mer sequences for MHC Class I epitope prediction",
}



sub execute {
    my $self = shift;
    #my $tmp_dir = Genome::Sys->create_temp_directory();

    my $input_fh = Genome::Sys->open_file_for_reading($self->input_file);
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    
    while (my $line = $input_fh->getline) {
		chomp $line;
		$line =~ s/[*]$//g;
		my @prot_arr =  split(/\t/, $line);
		if ( $prot_arr[15] =~ /^p.([A-Z])(\d+)([A-Z])/ && $prot_arr[13] eq $self->trv_type)
		{
			open(OUT, '>>', $self->output_file) or die $!;
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
         		if ($position < 8)
        		{
           			@wt_arr = @arr_wt_seq[ 0 ... 16];
           			$arr_wt_seq[$position]=$mt_aa;
           			@mt_arr = @arr_wt_seq[ 0 ... 16];
		     		print OUT ">WT.".$prot_arr[6].".".$prot_arr[15]."\n";
           			print OUT ( join "", @wt_arr);
           			print OUT "\n";
           			print OUT ">MT.".$prot_arr[6].".".$prot_arr[15]."\n";
           			print OUT ( join "", @mt_arr);
           			print OUT "\n";
        		}  
        		elsif ($position > ($#arr_wt_seq -8))
        		{	
           			@wt_arr = @arr_wt_seq[ $#arr_wt_seq -17 ... $#arr_wt_seq];
           			$arr_wt_seq[$position]=$mt_aa;
           			@mt_arr = @arr_wt_seq[ $#arr_wt_seq -17 ... $#arr_wt_seq];
           			print OUT ">WT.".$prot_arr[6].".".$prot_arr[15]."\n";
           			print OUT ( join "", @wt_arr);
           			print OUT "\n";
           			print OUT ">MT.".$prot_arr[6].".".$prot_arr[15]."\n";
           			print OUT ( join "", @mt_arr);
           			print OUT "\n";
        		}	
       		 	elsif (($position >= 8) && ($position  <= ($#arr_wt_seq -8)))
        		{
           			@wt_arr = @arr_wt_seq[ $position-8 ... $position+8];
           			$arr_wt_seq[$position]=$mt_aa;
           			@mt_arr = @arr_wt_seq[ $position-8 ... $position+8];
           			print OUT ">WT.".$prot_arr[6].".".$prot_arr[15]."\n";
           			print OUT ( join "", @wt_arr);
           			print OUT "\n";
           			print OUT ">MT.".$prot_arr[6].".".$prot_arr[15]."\n";
           			print OUT ( join "", @mt_arr);
           			print OUT "\n";
       				}
        		else 
        		{
        		print OUT "NULL"."\t".$position."\n";
        		}
        		OUT->close;
        		
    }
}  
    }
1;
}
return 1;
