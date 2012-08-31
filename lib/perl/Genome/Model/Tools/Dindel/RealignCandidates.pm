package Genome::Model::Tools::Dindel::RealignCandidates;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel::RealignCandidates {
    is => 'Command',
    has => [
    variant_file=> {
        is=>'String',
        is_input=>1,
    },
    output_file=>{
        is=>'String',
        is_input=>1,
    },
    ref_fasta=> {
        is=>'String',
        is_input=>1,
    },
    ],
};

sub help_brief {
    'Left-shift indels you got from a vcf'
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}


sub execute {
    my $self = shift;
    my $dindel_location = "/gscmnt/gc2146/info/medseq/dindel/binaries/dindel-1.01-linux-64bit";
    my $ref = $self->ref_fasta;
    my $output = $self->output_file;
    my $input = $self->variant_file;
    my $cmd = "$dindel_location --analysis realignCandidates --varFile $input --outputFile $output --ref $ref";
    return Genome::Sys->shellcmd(cmd=>$cmd);
}

1;
