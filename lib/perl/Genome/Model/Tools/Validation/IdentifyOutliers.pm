package Genome::Model::Tools::Validation::IdentifyOutliers;

use strict;
use warnings;

use IO::File;
use IO::Handle;
use Genome;
use List::MoreUtils

my %stats = ();

class Genome::Model::Tools::Validation::IdentifyOutliers {
    is => 'Command',

    has => [
        input_file => { 
            is => 'Text', 
            doc => "annotation format with readcounts added and headers required. Often named snvs.indels.annotated",
            is_optional => 0,
        },
        output_file => { 
            is => 'Text', 
            doc => "file into which high-confidence output should be placed",
            is_optional => 0,
        },
        lq_output_file => { 
            is => 'Text', 
            doc => "file into which low-confidence output should be placed",            
            is_optional => 1,
        },
        error_rate => { 
            is => 'Number', 
            doc => "The assumed per site error rate. Should be between 0 and 1", 
            is_optional => 1, 
            default => '0.01'
        },
        tumor_purity => { 
            is => 'Number', 
            doc => "The estimated purity of the tumor. Should be between 0 and 1", 
            is_optional => 1, 
            default => '1.0',
            is_optional => 1,
        },
        normal_contamination_rate => { 
            is => 'Number', 
            doc => "The estimated contamination rate of the normal with the tumor. Should be between 0 and 1", is_optional => 1, 
            default => '0.0',
            is_optional => 1,
        },
        llr_cutoff => { 
            is => 'Number', 
            doc => "minimum log likelihood to report a potential miscall",
            default => 10,
            is_optional => 1,
        },
        # min_p_value => { 
        #     is => 'Number', 
        #     doc => "minimum somatic confidence of Varscan to report the potential miscall", 
        #     default => 0.001
        # },
        # print_all => {
        #     is => 'Boolean', 
        #     doc => 'print information on all sites regardless of cutoffs', 
        #     default => 0,
        #     is_optional => 1,
        # },
        
        print_only_somatic => {
            is => 'Boolean', 
            doc => 'only print somatic sites, remove non-somatic', 
            default => 0,
            is_optional => 1,
        }

        ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Identify potential outliers in validation data for manual review";
}

sub help_synopsis {
    return <<EOS
This script attempts to identify SNV sites that are outliers with respect to the known tumor and normal contamination rates. Sites of low confidence or contradictory calls are printed to standard out.
EOS
}

sub execute {
    my $self = shift;
    my $fh = IO::File->new($self->input_file,"r");
    unless($fh) {
        $self->error_message("Unable to open input");
        return;
    }

    open(OUTFILE, ">" . $self->output_file);
    if(defined($self->lq_output_file)){
        open(LQFILE, ">" . $self->lq_output_file);
    }

    my $normal_ref_col;
    my $normal_var_col;
    my $tumor_ref_col;
    my $tumor_var_col;

    #scan through file and identify where we think there are outliers
    while(my $line = $fh->getline) {
        chomp $line;
        my @F = split /\t/, $line;
        
        #get fields to query from header
        if($line =~ /chromosome_name/){
            $normal_ref_col = List::MoreUtils::first_index {$_ eq "Normal_ref_count" } @F;
            $normal_var_col = List::MoreUtils::first_index {$_ eq "Normal_var_count" } @F;
            $tumor_ref_col = List::MoreUtils::first_index {$_ eq "Tumor_ref_count" } @F;
            $tumor_var_col = List::MoreUtils::first_index {$_ eq "Tumor_var_count" } @F;

            print OUTFILE $line . "\tFilterCall\tProb\tLLR\n";
            if(defined($self->lq_output_file)){
                print LQFILE $line . "FilterCall\tProb\tLLR\n";
            }
            next;
        }

        
        my $normal_ref = $F[$normal_ref_col];
        my $normal_var = $F[$normal_var_col];
        my $tumor_ref = $F[$tumor_ref_col];
        my $tumor_var = $F[$tumor_var_col];


        #my ($normal_ref, $normal_var, $tumor_ref,$tumor_var,$varscan_pvalue,$varscan_call) = @fields[4,5,8,9,14,12];
#        next if($varscan_pvalue > $self->min_p_value);

        my $llr = 0;
        my $max_llr;
        my $max_call;
        my $max2_llr;
        my $max2_call;
        my $marginal_probability;
        
        #pass through long indels because we can't grab readcounts.
        if($normal_ref eq "NA" || $normal_var eq "NA" || $tumor_ref eq "NA" || $tumor_var eq "NA"){
            print OUTFILE $line . "\tSomatic\tNA\tNA\n";
            next;
        }
        

        if($normal_ref + $normal_var > 0 && $tumor_ref + $tumor_var > 0) {

            #calculate LLR
            #previous attempt was a miserable failure due to degrees of freedom
            #use binomial model where variant supporting read is success!
            #take error rate as 1/1000 just for kicks and so we can actually calculate the results. Other numbers may be more appropriate
            #Want to test several models and pick the most likely one
            #Reference: normal expectation is 0.001 and tumor expectation is 0.01
            #Germline Het: normal expectation is 0.5 and tumor expectation is 0.5
            #Somatic Het: normal expectation is 0.001 and tumor expectation in 0.5
            #Germline Homozygote: normal expectation is 0.999 and tumor expectation is 0.999
            #LOH (variant): germline is 0.5 and tumor is 0.999
            #LOH (ref): germline is 0.5 and tumor is 0.001
            my $error_rate = $self->error_rate;

            my $normal_contamination_rate = $self->normal_contamination_rate;
            my $tumor_purity = $self->tumor_purity;

            my $homozygous_expect = 1 - $error_rate;
            my $heterozygous_expect = 0.5;


            #weighted average of the frequencies
            my $error_expect = ($tumor_var + $normal_var)/($tumor_var + $tumor_ref + $normal_var + $normal_ref); 
            #my $error_expect = ($tumor_var/($tumor_var + $tumor_ref) + $normal_var/($normal_var + $normal_ref))/2; #average of the frequencies


            $error_expect ||= $error_rate;
            if($error_expect == 1) {
                $error_expect = 1 - $error_rate;
            }

            my $tumor_freq = $tumor_var/($tumor_var + $tumor_ref);
            if($tumor_freq == 0) {
                $tumor_freq = $error_rate;
            }
            elsif($tumor_freq == 1) {
                $tumor_freq -= $error_rate;
            }
            my %calls = (   
                Reference       => [$error_rate, $error_rate],
                Germline_het    => [$heterozygous_expect, $heterozygous_expect],
                Germline_hom    => [$homozygous_expect, $homozygous_expect],
                Somatic     => [$error_rate + ($tumor_freq / $tumor_purity * $normal_contamination_rate), $tumor_freq],
                LOH_variant     => [$heterozygous_expect, $homozygous_expect],
                LOH_ref         => [$heterozygous_expect, $error_rate],
                NotSomatic           => [$error_expect, $error_expect], #In general this fits too well. It would need some other parameter that competed equally well to work. Perhaps, for instance, the existing tumor frequency and error in the normal 
            );

            my $llr_calculator = $self->generate_llr_calculator($normal_ref,$normal_var,$tumor_ref,$tumor_var);

            for my $call (keys %calls) { 
                my $llr = $llr_calculator->(@{$calls{$call}});
                #add to marginal
                unless(defined $marginal_probability) {
                    $marginal_probability = $llr;
                }
                else {
                    if($llr >= $marginal_probability) {
                        $marginal_probability = $llr + log(1 + exp($marginal_probability-$llr));
                    }
                    else {
                        $marginal_probability = $marginal_probability + log(1 + exp($llr-$marginal_probability));
                    }
                }

                unless(defined $max_llr && $llr < $max_llr) {
                    $max2_call = $max_call;
                    $max2_llr = $max_llr;
                    $max_llr = $llr;
                    $max_call = $call;

                }
                else {
                    unless(defined $max2_llr && $llr < $max2_llr) {
                        $max2_call = $call;
                        $max2_llr = $llr;
                    }
                }
            }
            my $llr = defined $max_llr ? $max_llr-$max2_llr : 0;

            # if((($varscan_call eq 'Somatic' && $max_call ne 'Somatic') || 
            #     (($varscan_call ne 'Somatic' && $varscan_call ne 'LOH')  && $max_call eq 'Somatic')) || 
            #    ($llr < $self->llr_cutoff && ($varscan_call eq 'Somatic')) || 
            #    $self->print_all) {
            my $call = (defined $max_call ? $max_call : '-');

            if(($llr < $self->llr_cutoff)){
                if(defined($self->lq_output_file)){
                    if(($self->print_only_somatic && ($call eq "Somatic")) || !($self->print_only_somatic)){
                        print LQFILE $line;
                        print LQFILE "\t" . $call;
                        print LQFILE "\t" . (defined $max_llr ? exp($max_llr-$marginal_probability) : '-');
                        print LQFILE "\t" . (defined $max_llr && defined $max2_llr ? $max_llr-$max2_llr : '-');
                        print LQFILE "\n";
                    }
                }
            } else {
                if( ($self->print_only_somatic && ($call eq "Somatic")) || !($self->print_only_somatic)){
                    print OUTFILE $line;
                    print OUTFILE "\t" . $call;
                    print OUTFILE "\t" . (defined $max_llr ? exp($max_llr-$marginal_probability) : '-');
                    print OUTFILE "\t" . (defined $max_llr && defined $max2_llr ? $max_llr-$max2_llr : '-');
                    print OUTFILE "\n";
                }
            }
            
        } else {
            # if there are no reads in either the tumor or normal, we set it to Non-Somatic with low probability
            # (no confidence if no reads)
            if (defined($self->lq_output_file)){
                unless($self->print_only_somatic){
                    print LQFILE $line;
                    print LQFILE "\tNotSomatic\t1\t0\n";
                }
            }
        }
        

    }
    close(OUTFILE);
    if(defined($self->lq_output_file)){
        close(LQFILE);
    }
    return 1;

}

sub generate_llr_calculator {
    my ($self,$normal_ref_reads,$normal_var_reads,$tumor_ref_reads,$tumor_var_reads) = @_;
    return sub {
        my ($normal_expect, $tumor_expect) = @_;
        my $not_normal_expect = 1 - $normal_expect;
        my $not_tumor_expect = 1 - $tumor_expect;
        return ($normal_ref_reads * log($not_normal_expect) + $normal_var_reads * log($normal_expect) + $tumor_ref_reads*log($not_tumor_expect) + $tumor_var_reads*log($tumor_expect));
    }
}

1;

