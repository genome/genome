package Genome::Model::Tools::Dindel::VcfToDindel;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel::VcfToDindel {
    is => 'Command',
    has => [
    input_vcf=> {
        is=>'String',
        is_input=>1,
    },
    output_dindel_file=>{
        is=>'String',
        is_input=>1,
        is_output=>1,
    },
    ref_fasta=> {
        is=>'String',
        is_input=>1,
    },
    ],
};

sub help_brief {
    'Turn vcf into stupid dindel format'
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
    my $script_location = "/gscmnt/gc2146/info/medseq/dindel/dindel-1.01-python/convertVCFToDindel.py";
    my $ref = $self->ref_fasta;
    my $output = $self->output_dindel_file;
    my $input = $self->input_vcf;
    my $cmd = "python $script_location --inputFile $input --outputFile $output --refFile $ref";
    return Genome::Sys->shellcmd(cmd=>$cmd);
}

1;
