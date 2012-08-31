package Genome::Model::Tools::Analysis::Indels::FilterSomatic;

use strict;
use warnings;

use IO::File;
use IO::Handle;
use Genome;

class Genome::Model::Tools::Analysis::Indels::FilterSomatic {
    is => 'Command',

    has => [
    somatic_file    => { is => 'Text', doc => "Indel results from gmt analysis indels compile-somatic" },
    fdr => { is => 'Number', doc => "The maximum fdr value to pass the filter", default => 0.000001 },
    llr => { is => 'Number', doc => "The minimum llr difference to pass filter", is_optional => 1, default => '10'},
    has_error_column => { is => 'Boolean', doc => 'Whether or not the annotation file has an error column', default => 1},
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
    my $error_col_adjustment = $self->has_error_column ? 0 : 1;

    #check the inputs before doing any significant work
    my $fh = IO::File->new($self->somatic_file,"r");
    unless($fh) {
        $self->error_message("Unable to open ". $self->somatic_file);
        return;
    }

    while(my $line = $fh->getline) {
        chomp $line;
        next if $line =~ /call_llr_to_second_best/;

        my @fields = split /\t/, $line;
        next if($fields[33 - $error_col_adjustment ] eq '-' || $fields[31 - $error_col_adjustment] eq '-');
        next if($fields[33 - $error_col_adjustment] > $self->fdr);
        next if($fields[31 - $error_col_adjustment] < $self->llr);
        print $line,"\n";
    }

    return 1;

}

1;

