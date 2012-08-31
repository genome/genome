package Genome::Model::Tools::Snp::Filters::DtrSeeFourFive;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Workflow;
use FileHandle;

class Genome::Model::Tools::Snp::Filters::DtrSeeFourFive{
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
        c_file =>
        {
            is => 'String',
            is_input => 1,
            doc => 'c4.5 tree output'
        },
        debug_mode => {default => 0},
        probability_mode => {default => 0},
        data_format =>
        {
            type => 'String',
            is_optional => 0,
            doc => 'Which Maq::Metrics::Dtr module to use',
            default => 'MaqOSixThree',

        },

        ]
};

sub help_synopsis {
    "Tool version of Dtr3e"
}

sub execute {
    my $self=shift;
    my $file = $self->experimental_metric_model_file;

    my $type = $self->data_format;
    my $dtr = eval "Genome::Model::Tools::Maq::Metrics::Dtr::$type->create()";
    unless(defined($dtr)) {
        $self->error_message($@);
        return;
    }
#    my $dtr_file = Genome::Model::Tools::Maq::Metrics::ConvertForDtr->execute(input=>$self->experimental_metric_model_file)->result;
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
    $handle->open($file, "r") or die "Couldn't open experimental metrics file\n";
    #my $dtr_handle = new FileHandle;
    #$dtr_handle->open($dtr_file, "r") or die "Couldn't open indendified metrics file\n";
    
    #we want to pass over dtr's columns but save the other file.
    #my $header_line = $dtr_handle->getline; #store header
    my $other_head= $handle->getline;
    chomp $other_head;
    #chomp $header_line;
    my @headers = $dtr->headers; #split(/,\s*/,$header_line);
    for (@headers) { s/\-/MINUS/g; s/\+/PLUS/g; s/\s/SPACE/; };

    
    my $keep_handle = new FileHandle;
    my $remove_handle = new FileHandle;
    my $keep_file = $self->decision_tree_output;
    my $remove_file = $self->basedir . 'chr' . $self->ref_seq_id . '.remove.csv';
    $keep_handle->open("$keep_file","w") or die "Couldn't open keep output file $keep_file: $!\n";
    $remove_handle->open("$remove_file","w") or die "Couldn't open remove output file $remove_file: $!\n";

    print $keep_handle $other_head ,"\n";
    print $remove_handle $other_head , "\n";

    #we could probably take the below method calls and turn them into an execute!!!! OMG
    my $c45_object = Genome::Model::Tools::SeeFourFive::Tree ->create();
    if($self->debug_mode) {
        $c45_object->debug_mode(1);
    }
    if($self->probability_mode) {
        $c45_object->probability_mode(1);
    }
    $c45_object->c45_file($self->c_file);
    $c45_object->load_trees;
  
    my $fref = $c45_object->generate_callback_for_headers(@headers);
    unless ($fref) {
        $self->error_message("No classifier from c5 trial");
        return;
    }
    
    while(my $line=$handle->getline) {
        chomp $line;
        
        next if ($line =~ /^chromosome/);
        #next if $line =~ /N/; #This was causing a bug with the snpfilter output
        
        my %data; 
        my @dtr_data = $dtr->make_attribute_array($line);
        unless(@dtr_data) {
            next;
        }
        @data{@headers} = @dtr_data;
        my ($decision, $leaf_number);
        my @answers;
        if($self->debug_mode || $self->probability_mode) {
            @answers =  $fref->(\%data);
            $decision = shift @answers;
        }
        else {
            ($decision) = $fref->(\%data);
        }

        unless ($decision) {
            ($decision) = $fref->(\%data);
            die "failed to get a decision for $line!";
        }
       
        ##this set of ifs is pretty ugly...refactor when we care enough
        if ($decision eq 'G') {
            if($self->debug_mode || $self->probability_mode) {
                #my $dtr_line = join ",",@data{@headers};
                print $keep_handle $line, "," , join(',',@answers) ,"\n";
            }
            else {
                print $keep_handle $line, "\n";
            }
        }  
        elsif ($decision eq 'WT') {
             if($self->debug_mode || $self->probability_mode) {
                my $dtr_line = join ",",@data{@headers};
                print $remove_handle $line, ",", join(',',@answers) ,"\n";
            }
            else {
                print $remove_handle $line, "\n";
            }
        }
        else {
            die "unhandled decision $decision from line $line!";
        }
    }
    $keep_handle->close();
    $remove_handle->close();
    $handle->close();
    return 1;
}

sub _demo_c5src {
    return <<EOS

SeeFive.0 [Release 2.05]     Wed Aug  6 18:16:05 2008
-------------------

    Options:
        Application `./testing/first'
        Rule-based classifiers
        Boosted classifiers

Read 870 cases (21 attributes) from ./testing/first.data
Read misclassification costs from ./testing/first.costs

-----  Trial 0:  -----

Rules:

Rule 0/1: (279/24, lift 2.1)
        avg_sum_of_mismatches <= 35
        avg_base_quality > 18
        unique_variant_for_reads_ratio > 0.08
        ->  class G  [0.911]

Rule 0/2: (216/29, lift 2.0)
        pre-27_unique_variant_context_reads_ratio > 0.3333333
        ->  class G  [0.862]

Rule 0/3: (187/5, lift 1.7)
        max_base_quality <= 26
        ->  class WT  [0.968]

Rule 0/4: (220/10, lift 1.7)
        avg_base_quality <= 16
        ->  class WT  [0.950]

Rule 0/5: (98/5, lift 1.7)
        maq_snp_quality <= 23
        ->  class WT  [0.940]

Rule 0/6: (101/25, lift 1.3)
        max_mapping_quality <= 45
        ->  class WT  [0.748]

Rule 0/7: (654/191, lift 1.3)
        pre-27_unique_variant_context_reads_ratio <= 0.3333333
        ->  class WT  [0.707]

Default class: WT

-----  Trial 1:  -----

Rules:

Rule 1/1: (406.9/25.2, lift 1.6)
        maq_snp_quality > 35
        avg_mapping_quality > 31
        max_base_quality > 26
        unique_variant_for_reads_ratio > 0.02222222
        unique_variant_rev_reads_ratio > 0.05263158
        ->  class G  [0.936]



EOS
}

=cut

        my ($chr,$position,$al1,$al2,$qvalue,$al2_read_hg,$avg_map_quality,$max_map_quality,$n_max_map_quality,$avg_sum_of_mismatches,$max_sum_of_mismatches,$n_max_sum_of_mismatches, $base_quality, $max_base_quality, $n_max_base_quality, $avg_windowed_quality, $max_windowed_quality, $n_max_windowed_quality, $for_strand_unique_by_start_site, $rev_strand_unique_by_start_site, $al2_read_unique_dna_context, $for_strand_unique_by_start_site_pre27,$rev_strand_unique_by_start_site_pre27,$al2_read_unique_dna_context_pre27) = split ", ", $line;
        my $al2_read_unique_dna_start = $for_strand_unique_by_start_site + $rev_strand_unique_by_start_site;
        my $al2_read_unique_dna_start_pre27 = $for_strand_unique_by_start_site_pre27 + $rev_strand_unique_by_start_site_pre27;

=cut

1;
