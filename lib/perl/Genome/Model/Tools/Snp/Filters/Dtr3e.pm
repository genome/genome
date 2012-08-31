package Genome::Model::Tools::Snp::Filters::Dtr3e;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Workflow;

class Genome::Model::Tools::Snp::Filters::Dtr3e {
    is => 'Command',
    has => [
    experimental_metric_model_file => 
    {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'File of experimental metrics for the model'
    },
    basedir => 
    { 
        type => 'String',
        is_optional =>0,
        is_input => 1,
        doc => 'Place and prefix for files'
    },
    specificity =>
    {
        type => 'String',
        is_optional => 1,
        doc => 'its the specificity' , 
        default=> 'default'
    },
    ref_seq_id => 
    { 
        is => 'String',
        is_input => 1,
        doc => 'Put some docmumentation here' 
    },
    decision_tree_output => 
    {
        doc => ".keep output file for this step.", 
        is => 'String',
        is_output => 1,
        calculate => q| 
             return $self->basedir . 'chr' . $self->ref_seq_id . '.keep.csv';
        |
    },
    ]
};

sub help_synopsis {
    "Tool version of Dtr3e"
}

sub execute {
    my $self=shift;
    my $file = $self->experimental_metric_model_file;
    my $basename=$self->basedir . "/filtered";
    my $specificity = $self->specificity; 
    my %specificity_maqq = (
        min => [ 5, 16 ],
        default => [ 30, 16 ],
        90 => [ 40, 22 ],
        95 => [ 50, 26 ],
        max => [ 120, 26 ],
    );


    my $spec_ref = (exists($specificity_maqq{$specificity})) ?
    $specificity_maqq{$specificity} : $specificity_maqq{default};
    my ($qvalue_level, $bq) = @$spec_ref;
    $qvalue_level ||= 30;
    $bq ||= 16;




    my $handle = new FileHandle;
    $handle->open($file, "r") or die "Couldn't open annotation file\n";

    my $header_line = $handle->getline; #store header
    my $keep_handle = new FileHandle;
    my $remove_handle = new FileHandle;
    my $keep_file = $self->decision_tree_output;
    my $remove_file = $basename . 'chr' . $self->ref_seq_id . '.remove.csv';
    $keep_handle->open("$keep_file","w") or die "Couldn't open keep output file\n";
    $remove_handle->open("$remove_file","w") or die "Couldn't open remove output file\n";

    chomp $header_line;
    print $keep_handle $header_line . ", rule\n";
    print $remove_handle $header_line . ", rule\n";


    my %result;
#print new header
while(my $line=$handle->getline) {
    chomp $line;
    $line =~ s/NULL/0/g;
    next if $line =~ /^chromosome/;    
    my ($chr,$position,$al1,$al2,$qvalue,$al2_read_hg,$avg_map_quality,$max_map_quality,$n_max_map_quality,$avg_sum_of_mismatches,$max_sum_of_mismatches,$n_max_sum_of_mismatches, $base_quality, $max_base_quality, $n_max_base_quality, $avg_windowed_quality, $max_windowed_quality, $n_max_windowed_quality, $for_strand_unique_by_start_site, $rev_strand_unique_by_start_site, $al2_read_unique_dna_context, $for_strand_unique_by_start_site_pre27,$rev_strand_unique_by_start_site_pre27,$al2_read_unique_dna_context_pre27) = split ", ", $line;
    my $al2_read_unique_dna_start = $for_strand_unique_by_start_site + $rev_strand_unique_by_start_site;
    my $al2_read_unique_dna_start_pre27 = $for_strand_unique_by_start_site_pre27 + $rev_strand_unique_by_start_site_pre27;

    next if $qvalue < 15; #hardcoded to make this work, filtering out anything that is not q 15

    my $decision = 'keep';
    my $rule = 'none';

    if ($al2_read_unique_dna_start > 7 &&
    $qvalue >= $qvalue_level) {
        #Rule 5:
        #        # of unique genomic reads supporting variant allele(starting point) > 7
        #    ->  class G  [96.8%]
        #
        $decision = 'keep';
        $rule = '5';
#        } elsif ($al2_read_unique_dna_start > 2 &&
#                         $base_quality <= $bq + 2) {
#            #Rule 2:
#            #        # of unique genomic reads supporting variant allele(starting point) > 2
#            #        Base Quality > 18
#            #    ->  class G  [90.7%]
#            #
#            $decision = 'keep';
#            $rule = '2';
        } elsif ($max_base_quality <= 26) {
            #Rule 11:
            #        Max Base Quality <= 26
            #    ->  class WT  [85.2%]
            #
            $decision = 'remove';
            $rule = '11';
        } elsif ($al2_read_unique_dna_start <= 7 &&
        $qvalue < $qvalue_level) {
            #Rule 7:
            #        # of unique genomic reads supporting variant allele(starting point) <= 7
            #        Maq SNP q-value <= 29
            #    ->  class WT  [81.5%]
            #
            $decision = 'remove';
            $rule = '7';
        } elsif ($al2_read_unique_dna_start <= 2) {
            #Rule 8:
            #        # of unique genomic reads supporting variant allele(starting point) <= 2
            #    ->  class WT  [74.7%]
            #
            $decision = 'remove';
            $rule = '8';
        } elsif ($base_quality <= $bq) {
            #Rule 3:
            #        Base Quality <= 16
            #    ->  class WT  [73.0%]
            #
            $decision = 'remove';
            $rule = '3';
        } else {
            #Default class: G
            $decision = 'keep';
            $rule = 'default';
        }

        ###additional arbitrary filter on pre-27 readcounts
        #unless($al2_read_unique_dna_start_pre27 > 2) {
        #    $decision = 'remove';
        #    $rule = 'pre-27';
        #}

        if ($decision eq 'keep') {
            print $keep_handle $line,", $rule\n";    
        }  else {
            print $remove_handle $line,", $rule\n";    
        }
    }
    $keep_handle->close();
    $remove_handle->close();
    $handle->close();
    return 1;
}



1;
