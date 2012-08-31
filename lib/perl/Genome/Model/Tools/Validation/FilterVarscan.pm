package Genome::Model::Tools::Validation::FilterVarscan;

use strict;
use warnings;

use IO::File;
use IO::Handle;
use Genome;

class Genome::Model::Tools::Validation::FilterVarscan {
    is => 'Command',

    has => [
    varscan_file    => { is => 'Text', doc => "varscan results" },
    min_coverage => { is => 'Number', doc => "minimum number of reads total for both tumor and normal", default => 30 },
    max_p_value => { is => 'Number', doc => "minimum somatic confidence of Varscan", default => 0.001},
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Filter to confident varscan calls";
}

sub help_synopsis {
    return <<EOS
EOS
}

sub execute {
    my $self = shift;
    my $fh = IO::File->new($self->varscan_file,"r");
    unless($fh) {
        $self->error_message("Unable to open varscan file");
        return;
    }

    while(my $line = $fh->getline) {
        chomp $line;
        my @fields = split /\t/, $line;

        my ($normal_ref, $normal_var, $tumor_ref,$tumor_var,$varscan_pvalue,$varscan_call) = @fields[4,5,8,9,14,12];
        next if($varscan_pvalue > $self->max_p_value);
        next if(($normal_ref + $normal_var) < $self->min_coverage || ($tumor_ref + $tumor_var) < $self->min_coverage);
        print $line, "\n";
    }

    return 1;

}

1;

