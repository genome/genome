package Genome::Model::Tools::Dindel::AnalyzeWindowFile;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::Dindel::AnalyzeWindowFile {
    is => 'Command',
    has => [
    window_file=> {
        is=>'String',
        is_input=>1,
    },
    library_metrics_file=>{
        is=>'String',
        is_input=>1,
        doc=>'from step one... getCigarIndels',
    },
    bam_file=> {
        is=>'String',
        is_input=>1,
    },
    output_prefix=> {
        is=>'String',
        is_input=>1,
        is_output=>1,
    },
    ref_fasta=> {
        is=>'String',
        is_input=>1,
    },
    output_bam=> {
        is=>'Number',
        default=>0,
        is_optional=>1,
        is_input=>1,
    },
    ],
};

sub help_brief {
    'Actual slow part of dindel-- analysis'
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
    my (undef, $callback_script_path) = File::Basename::fileparse(__FILE__);
    my $callback_script = $callback_script_path . "merge_bam_callback.pl";

    if($self->output_bam) {
        unless(-s $callback_script) {
            $self->error_message("unable to locate dindel helper script at location: $callback_script\n");
            return undef;
        }
    }
    my $ref = $self->ref_fasta;
    my $output = $self->output_prefix;
    my $input = $self->window_file;
    my $lib_file = $self->library_metrics_file;
    my $bam = $self->bam_file;
    my $cmd = "$dindel_location --analysis indels --doDiploid --bamFile $bam --varFile $input --outputFile $output --ref $ref --libFile $lib_file";
    if($self->output_bam) {
        $cmd .= " --outputRealignedBAM";
        $cmd .= " --processRealignedBAM $callback_script";
    }
    return Genome::Sys->shellcmd(cmd=>$cmd);
}

1;
