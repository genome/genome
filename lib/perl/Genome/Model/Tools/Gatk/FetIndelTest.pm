package Genome::Model::Tools::Gatk::FetIndelTest;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Spreadsheet::WriteExcel;
use Sort::Naturally qw( nsort );
use Genome::Model::Tools::Varscan;

class Genome::Model::Tools::Gatk::FetIndelTest {
    is => 'Command',
    has => [
    gatk_verbose_output => { 
        type => 'String',
        is_optional => 0,
        doc => "the gatk output to test",
    },
    p_value_cutoff => { 
        type => 'Float',
        is_optional => 0,
        doc => "the maximum p-value to report a line",
        default => '0.001',
    },
    ]
};



sub execute {
    my $self=shift;

    my $indel_call_file = $self->gatk_verbose_output;

    my $fh = IO::File->new($indel_call_file,"r");
    unless($fh) {
        $self->error_message("Unable to open $indel_call_file: $!");
        return;
    }
    while(my $line = $fh->getline) {
        chomp $line;
        next if $line =~ /REDUCE RESULT/;
        my @fields = split /\t/, $line;
        my ($normal_counts_field) = grep(/^N_OBS_COUNTS/,@fields);
        my ($tumor_counts_field) = grep(/^T_OBS_COUNTS/,@fields);

        my ($normal_var,$normal_total) = $normal_counts_field =~ /:(\d+)\/\d+\/(\d+)/;
        my $normal_ref = $normal_total - $normal_var;   #this would not count other indel alleles 
        my ($tumor_var,$tumor_total) = $tumor_counts_field =~ /:(\d+)\/\d+\/(\d+)/;
        my $tumor_ref = $tumor_total - $tumor_var;   #this would not count other indel alleles 

        my $p_value = Varscan::FisherTest::calculate_p_value($normal_ref, $normal_var, $tumor_ref, $tumor_var);
        print $line,"\t",$p_value,"\n" if($self->p_value_cutoff >= $p_value);
    }        

    return 1;
}


1;

sub help_brief {
    "Generates tier1-4 indels"
}

sub help_detail {
    <<'HELP';
Hopefully this script will run the ucsc annotator on indels and then tier them for an already completed somatic model. Since this is done in the data directory, the model is then reallocated.
HELP
}

