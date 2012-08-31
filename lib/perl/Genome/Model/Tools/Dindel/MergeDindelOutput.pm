package Genome::Model::Tools::Dindel::MergeDindelOutput;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel::MergeDindelOutput {
    is => 'Command',
    has => [
    dindel_file_output_list=>{
        is=>'String',
        is_input=>1,
        doc=>'this should be a file of result filenames that dindel output in analyze-window-file, one filename/path per line',
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
    'make a vcf from dindel result files',
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
    my $script_location = "/gscmnt/gc2146/info/medseq/dindel/dindel-1.01-python/mergeOutputDiploid.py";
    my $output = $self->output_file;
    my $input = $self->dindel_file_output_list;
    my $ref = $self->ref_fasta;
    my $cmd = "python $script_location --inputFiles $input --outputFile $output --refFile $ref";
    return Genome::Sys->shellcmd(cmd=>$cmd);
}

1;
