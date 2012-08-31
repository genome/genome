package Genome::Model::Tools::Dindel::MakeDindelWindows;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel::MakeDindelWindows {
    is => 'Command',
    has => [
    input_dindel_file=>{
        is=>'String',
        is_input=>1,
        doc=>'file of dindel formatted indels to examine. get this from getcigarindels or vcftodindel followed by realigncandidates',
    },
    output_prefix=>{
        is=>'String',
        is_input=>1,
    },
    num_windows_per_file=> {
        is=>'String',
        is_input=>1,
        is_optional=>1,
        default=>1000,
    },
    ],
};

sub help_brief {
    'make window files for dindel parallelization'
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
    my $script_location = "/gscmnt/gc2146/info/medseq/dindel/dindel-1.01-python/makeWindows.py";
    my $output = $self->output_prefix;
    my $input = $self->input_dindel_file;
    my $num_windows = $self->num_windows_per_file;
    my $cmd = "python $script_location --inputVarFile $input --windowFilePrefix $output --numWindowsPerFile $num_windows";
    return Genome::Sys->shellcmd(cmd=>$cmd);
}

1;
