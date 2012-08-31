package Genome::Model::Tools::Analysis::Indels::CompileSomatic;

#This is mostly ported from CompileBowtieResults

use strict;
use warnings;

use IO::File;
use IO::Handle;
use Genome;
use Genome::Model::Tools::Varscan;
my %stats = ();

class Genome::Model::Tools::Analysis::Indels::CompileSomatic {
    is => 'Command',

    has => [
    tumor_file    => { is => 'Text', doc => "Indel results from gmt analysis indels compile-contig-counts" },
    normal_file => { is => 'Text', doc => "Indel results from gmt analysis indels compile-contig-counts" },
    output_file     => { is => 'Text', doc => "Output of indels with FET result", is_optional => 1 },
    error_rate => { is => 'Number', doc => "The assumed per site error rate. Should be between 0 and 1", is_optional => 1, default => '0.01'},
    tumor_purity => { is => 'Number', doc => "The estimated purity of the tumor. Should be between 0 and 1", is_optional => 1, default => '1.0'},
    normal_contamination_rate => { is => 'Number', doc => "The estimated contamination rate of the normal with the tumor. Should be between 0 and 1", is_optional => 1, default => '0.0'},
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Compile the tumor and normal results into one file and calculate a p-value";
}

sub help_synopsis {
    return <<EOS
This command compiles indel counts and calculates a p-value 
EXAMPLE:	gmt analysis indels compile-somatic --tumor-file [tumor.counts.tsv] --normal-file [normal.counts.tsv]
EOS
}

sub execute {
    my $self = shift;

    #check the inputs before doing any significant work
    my $tumor_fh = IO::File->new($self->tumor_file,"r");
    unless($tumor_fh) {
        $self->error_message("Unable to open ". $self->tumor_file);
        return;
    }

    my $normal_fh = IO::File->new($self->normal_file,"r");
    unless($normal_fh) {
        $self->error_message("Unable to open " . $self->normal_file);
        return;
    }

    my %indels; #hash to match up indels between files

    while(my $normal_line = $normal_fh->getline) {
        chomp $normal_line;

        my @fields = split /\t/, $normal_line;

        #generate unique identifier
        my $indel_name = uc(join("_",@fields[0..4]));

        if(exists($indels{$indel_name})) {
            $self->error_message("Multiple instances of the same indel present in file. Skipping line " . $normal_fh->input_line_number);
            next;
        }

        #store counts for later
        my ($ref, $var) = @fields[-3,-2];

        $indels{$indel_name} = {    line => $normal_line, 
            normal_ref => $ref,
            normal_var => $var,
        };
    }

    #print the header
    print join("\t",qw( chr stop start ref var type gene transcript species source version strand transcript_type mutation_type transcript_change amino_acid_change conservation domain all_domains deletion_annotation errors normal_total normal_ref normal_var normal_var_freq tumor_total tumor_ref tumor_var tumor_freq call call_prob call_llr_to_second_best FET_p_value FDR)),"\n";
    my %results;
    my $tested = 0;
    #read in tumor and calculate the p-value
    while(my $tumor_line = $tumor_fh->getline) {
        chomp $tumor_line;

        my @fields = split /\t/, $tumor_line;

        #generate unique identifier
        my $indel_name = uc(join("_",@fields[0..4]));

        unless(exists($indels{$indel_name})) {
            $self->error_message("Indel missing in the normal file or duplicated in the tumor file. Skipping line " . $tumor_fh->input_line_number);
            next;
        }

        #store counts for later
        my $entry = $indels{$indel_name};
        my ($ref, $var) = @fields[-3,-2];
        my $p_value = 1;

        my $llr = 0;
        my $max_llr;
        my $max_call;
        my $max2_llr;
        my $max2_call;
        my $marginal_probability;

        if(($entry->{normal_ref}+$entry->{normal_var}) > 0 && ($ref+$var) > 0) {
            $tested++;
            my $total_observed_reads = $ref + $var;
            my $expected_variant_reads = sprintf("%d",$entry->{normal_var}/($entry->{normal_ref} + $entry->{normal_var}) * $total_observed_reads);  #make sure it rounds properly to an integer
            my $expected_normal_reads = $total_observed_reads - $expected_variant_reads;

            $p_value = Varscan::FisherTest::calculate_p_value($expected_normal_reads, $expected_variant_reads, $ref, $var, 0);


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
            my $error_expect = ($var + $entry->{normal_var})/($total_observed_reads + $entry->{normal_ref} + $entry->{normal_var}); #weighted average of the two
            $error_expect ||= $error_rate;
            if($error_expect == 1) {
                $error_expect = 1 - $error_rate;
            }

            my $tumor_freq = $var/$total_observed_reads;
            if($tumor_freq == 0) {
                $tumor_freq = $error_rate;
            }
            elsif($tumor_freq == 1) {
                $tumor_freq -= $error_rate;
            }
            my %calls = (   Reference       => [$error_rate, $error_rate],
                Germline_het    => [$heterozygous_expect, $heterozygous_expect],
                Germline_hom    => [$homozygous_expect, $homozygous_expect],
                Somatic     => [$error_rate + ($tumor_freq / $tumor_purity * $normal_contamination_rate), $tumor_freq],
                LOH_variant     => [$heterozygous_expect, $homozygous_expect],
                LOH_ref         => [$heterozygous_expect, $error_rate],
                Error           => [$error_expect, $error_expect], #In general this fits too well. It would need some other parameter that competed equally well to work. Perhaps, for instance, the existing tumor frequency and error in the normal 
            );

            #my %calls = (   Reference       => [$error_rate, $error_rate],
            #    Germline_het    => [$heterozygous_expect, $heterozygous_expect],
            #    Germline_hom    => [$homozygous_expect, $homozygous_expect],
            #    Somatic_het     => [$error_rate, $heterozygous_expect],
            #    Somatic_hom     => [$error_rate, $homozygous_expect],
            #    LOH_variant     => [$heterozygous_expect, $homozygous_expect],
            #    LOH_ref         => [$heterozygous_expect, $error_rate],
            #    #Error           => [$error_expect, $error_expect], #In general this fits too well. It would need some other parameter that competed equally well to work. Perhaps, for instance, the existing tumor frequency and error in the normal 
            #);

            my $llr_calculator = $self->generate_llr_calculator($entry->{normal_ref},$entry->{normal_var},$ref,$var);

            for my $call (reverse sort keys %calls) {   #reverse sort ensures the Error status is only returned as a last resort...I hope
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
            $results{$indel_name} = { line => $entry->{line}."\t".join("\t",@fields[-4,-3,-2,-1],defined $max_call ? $max_call : '-',defined $max_llr ? exp($max_llr-$marginal_probability) : '-',defined $max_llr && defined $max2_llr ? $max_llr-$max2_llr : '-',), p => $p_value};

        }
        else {
            print STDOUT $entry->{line}."\t".join("\t",@fields[-4,-3,-2,-1],defined $max_call ? $max_call : '-',defined $max_llr ? exp($max_llr-$marginal_probability) : '-',defined $max_llr && defined $max2_llr ? $max_llr-$max2_llr : '-','-','-'),"\n";
        }

        delete $indels{$indel_name};
    }
    my $rank = 1;
    #calculate the FDR and finish it up
    for my $indel (sort { $results{$a}{p} <=> $results{$b}{p} } keys %results) {
        print $results{$indel}{line},"\t",$results{$indel}{p},"\t",$results{$indel}{p}*$tested/$rank++,"\n";
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

