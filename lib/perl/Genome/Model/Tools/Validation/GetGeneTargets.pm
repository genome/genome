package Genome::Model::Tools::Validation::GetGeneTargets;

use strict;
use warnings;
use FileHandle;
use Command;

class Genome::Model::Tools::Validation::GetGeneTargets {
    is => 'Command',
    has => [
    input_file => {
    type => 'String',
    is_optional => 0,
    doc => "input file containing gene list, one gene per line, without head",
    },
    output_file => {
    type => 'String',
    is_optional => 1,
    doc => "output file containing targets information, default is targets.tsv in your working directory",
    },   
    exon_file => {
    type => 'String',
    is_optional => 0,
    doc => "exon bed file with chromosome, start, stop, and gene name, separated by tab.  build37: /gscmnt/gc2108/info/medseq/ckandoth/bed_maker/NCBI-human.combined-annotation_58_37c_v2/all_CDS_and_ncRNA_24Chroms_Contigs_1BasedStart_2bpFlanks_MergedExons",
    },
    ] 
};

sub execute
{
        my $self=shift;
        $DB::single = 1;
        
        my $gene_file;
        my %genelist_hash;
        
        if ($self->input_file){
                $gene_file=$self->input_file;
                my $gene_fh=new FileHandle($gene_file);
                while(<$gene_fh>){
                chomp;
                my ($gene,@others)=split/\t/;
                $genelist_hash{$gene}=1;
                }
                
                my $summary_file;
                my $output_file;
                if ($self->output_file){
                        $output_file=$self->output_file;
                        $summary_file="$output_file.summary";
                }
                else{
                        $output_file="targets.tsv";
                        $summary_file="summary";
                }
                open OUT, ">$output_file";
                
                my $exon_file;
                if ($self->exon_file){
                        $exon_file=$self->exon_file;
                }
                my %exon_hash;
                my $exon_fh=new FileHandle($exon_file);
                while(<$exon_fh>){
                        chomp;
                        my $line=$_;
                        my ($chr,$start,$stop,$gene)=split/\t/;
                        my $pos=$chr."_".$start."_".$stop;
#                        if (defined $exon_hash{$gene}){  
                        if (defined $genelist_hash{$gene}){  
                                print OUT "$line\n";
                                push @{$exon_hash{$gene}}, $pos;
                        }
                        else{
                                my @arr=($pos);
                                $exon_hash{$gene}=\@arr;
                        }
               }
               close OUT;
               
               open SUM, ">$summary_file";
               my %selected_genes;

                print SUM "Summary: \n";
                print SUM "================\n";
                print SUM "Genes selected: \n";
                print SUM "================\n";
                print SUM "Gene\tnumber_of_exons\ttotal_bases\n";

                foreach my $gene (keys %genelist_hash){
                        if (defined $exon_hash{$gene}){
                                
                                my $total_bases=0;
                                my @exons=@{$exon_hash{$gene}};
                                my $num_exon=@exons;
                                foreach my $exon (@exons){
                                        my ($chr,$start,$stop)=split/\_/, $exon;
                                        my $bases=$stop-$start+1;
                                        print "please check this exon: $exon, the total bases <=0\n" if ($bases<=0);
                                        $total_bases+=$bases;
                                }
                                $selected_genes{bases}=$total_bases;
                                $selected_genes{exon_num}=$num_exon;
                                print SUM "$gene\t$num_exon\t$total_bases\n";
                                delete $genelist_hash{$gene};
                        }
                        
                }

                print SUM "\n\n==============================\n";
                print SUM "Warning: Genes without exons \n";
                print SUM "==============================\n";
                foreach my $g (keys %genelist_hash){
                        print SUM "$g\n";
                }


               close SUM;

        }
        else{
                print "Please enter an input file containing gene list, one gene per line\n";
                return 1;
        }
        
        
        return 1;
}
1;

sub help_brief {
    "gives out the exon regions for a given gene list"
}

sub help_detail {
    <<'HELP';
this script will print out the chr,start and end of exons for a given gene list, it will also generate a summary report containing the following: 
1) how many exons and bases per gene were selected
2) warn you about genes that had no exons in the file.
The summary file will be named your_output.summary along with your output_file, if output_file is not defined, then the summary file will be in your working directory, and named "summary". 
HELP
}
