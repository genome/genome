
package Genome::Model::Tools::Varscan::FilterVariantCalls;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Varscan::FilterVariantCalls    Process somatic pipeline output
#
#    AUTHOR:     Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#    CREATED:    12/09/2009 by D.K.
#    MODIFIED:   12/29/2009 by D.K.
#
#    NOTES:
#
#####################################################################################################################################

use strict;
use warnings;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Varscan::FilterVariantCalls {
    is => 'Command',

    has => [                                # specify the command's single-value properties (parameters) <---
        variant_file     => { is => 'Text', doc => "File containing varscan calls, e.g. status.varscan.snp" , is_optional => 0, is_input => 1},
        output_file  => { is => 'Number', doc => "P-value threshold for high confidence", is_optional => 1, is_input => 1, default_value => '0.07'},
        min_normal_freq => { is => 'Number', doc => "Minimum normal frequency", is_optional => 1, is_input => 1, default_value => '0'},
        max_normal_freq => { is => 'Number', doc => "Maximum normal frequency", is_optional => 1, is_input => 1, default_value => '100'},
        min_tumor_freq  => { is => 'Number', doc => "Minimum tumor frequency", is_optional => 1, is_input => 1, default_value => '20'},
        max_tumor_freq  => { is => 'Number', doc => "Maximum tumor frequency", is_optional => 1, is_input => 1, default_value => '100'},
        min_normal_cov  => { is => 'Number', doc => "Minimum normal coverage", is_optional => 1, is_input => 1, default_value => '10'},
        min_tumor_cov  => { is => 'Number', doc => "Minimum tumor coverage", is_optional => 1, is_input => 1, default_value => '10'},        
        skip_if_output_present => { is => 'Boolean', doc => "If set to 1, will not run if output is present" , is_optional => 1, is_input => 1, default_value => 0},
        output_file => { is => 'Text', doc => 'Output file for filtered calls', is_optional => 1, is_input => 1 },

    ],
    has_param => [
        lsf_resource => { value => 'select[tmp>1000] rusage[tmp=1000]'},
    ]
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Process output from Varscan somatic"                 
}

sub help_synopsis {
    return <<EOS
    This command processes output from Varscan somatic (.snp or.indel), classifying variants by somatic status
    (Germline/Somatic/LOH) and by confidence (high/low).
    EXAMPLE:    gmt varscan process-somatic ...
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;

    ## Get required parameters ##
    my $variant_file = $self->variant_file;
    my $min_normal_freq = $self->min_normal_freq;
    my $min_tumor_freq = $self->min_tumor_freq;
    my $max_normal_freq = $self->max_normal_freq;
    my $max_tumor_freq = $self->max_tumor_freq;
    my $min_normal_cov = $self->min_normal_cov;
    my $min_tumor_cov = $self->min_tumor_cov;
    my $output_file = $self->output_file;



    my %variants_by_status = ();
    $variants_by_status{'Somatic'} = $variants_by_status{'Germline'} = $variants_by_status{'LOH'} = '';

    open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

    ## Parse the variant file ##

    my $input = Genome::Sys->open_file_for_reading($variant_file);
    my $lineCounter = 0;

    while (<$input>) {
        chomp;
        my $line = $_;
        $lineCounter++;

        my @lineContents = split(/\t/, $line);

        if(($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name"))
        {
            #$file_header = $line;
        }
        else
        {
            my $offset = 0;
            ## If in formatted format (chrom start stop ref var), adjust the parsing offset ##
            my $test = $lineContents[2];
            $offset = 1 if($test =~ /[0-9]/);
            
            my $normal_cov = $lineContents[$offset + 4] + $lineContents[$offset + 5];
            my $normal_freq = $lineContents[$offset + 6];
            my $tumor_cov = $lineContents[$offset + 8] + $lineContents[$offset + 9];
            my $tumor_freq = $lineContents[$offset + 10];
            $normal_freq =~ s/\%//;
            $tumor_freq =~ s/\%//;           

            if($normal_cov >= $min_normal_cov && $tumor_cov >= $min_tumor_cov)
            {
                if($normal_freq >= $min_normal_freq && $normal_freq <= $max_normal_freq)
                {
                    if($tumor_freq >= $min_tumor_freq && $tumor_freq <= $max_tumor_freq)
                    {
                        print OUTFILE "$line\n";
                    }
                }
            }
        }
    }

    close($input);

    close(OUTFILE);

    return 1;
}


1;
