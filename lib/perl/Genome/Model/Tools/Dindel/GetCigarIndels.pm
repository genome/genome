package Genome::Model::Tools::Dindel::GetCigarIndels;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel::GetCigarIndels {
    is => 'Command',
    has => [
    input_bam=>{
        is=>'String',
        is_input=>1,
    },
    ref_fasta=>{
        is=>'String',
        is_input=>1,
    },
    output_prefix=>{
        is=>'String',
        is_input=>1,
    },
    ],
};

sub help_brief {
    'Run getCIGARindels'
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
    my $bam = $self->input_bam;
    my $ref = $self->ref_fasta;
    my $out_prefix = $self->output_prefix;
    my $cmd = "$dindel_location --analysis getCIGARindels --bamFile $bam --outputFile $out_prefix --ref $ref";
    return Genome::Sys->shellcmd(cmd=>$cmd);
}

1;
