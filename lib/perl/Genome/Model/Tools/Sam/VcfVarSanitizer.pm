package Genome::Model::Tools::Sam::VcfVarSanitizer;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;


class Genome::Model::Tools::Sam::VcfVarSanitizer {
    is  => 'Command',
    has => [
        var_file => {
            is  => 'String',
            doc => 'The input vcf var file',
        },
        out_file => {
            is  => 'String',
            doc => 'var output file after sanitizing, default is using the same var_file name',
            is_optional => 1,
        },
    ],
};


sub help_brief {
    'Sanitize samtools-mpileup var vcf output';
}

sub help_detail {
    return <<EOS
    Sanitize samtools-mpileup var output. Sometimes 4th column gets N as ref base and 5th column get . as alt
    
    16	35143303	.	N	T	19.8	.	DP=2;VDB=0.0190;AF1=1;AC1=2;DP4=0,0,1,1;MQ=46;FQ=-33	GT:PL:GQ	1/1:51,6,0:10
    16	44943302	.	N	T	45	.	DP=23;VDB=0.0000;AF1=1;AC1=2;DP4=0,0,23,0;MQ=14;FQ=-96	GT:PL:GQ	1/1:78,69,0:99
    17	22187134	.	N	G	29.5	.	DP=5;VDB=0.0374;AF1=1;AC1=2;DP4=0,0,0,4;MQ=29;FQ=-39	GT:PL:GQ	1/1:62,12,0:21

EOS
}


sub execute {
    my $self = shift;
    my $var_file = $self->var_file;

    unless (-s $var_file) {
        $self->error_message('Can not find valid SAM snp file: '.$var_file);
        return;
    }

    my $out_file;

    if ($self->out_file) {
        $out_file = $self->out_file;
    }
    else {
        my ($name, $path) = fileparse $var_file;
        $out_file = $path.$name.'.sanitize';
    }


    my $out_fh = Genome::Sys->open_file_for_writing($out_file) or return;
    my $var_fh = Genome::Sys->open_file_for_reading($var_file) or return;

    while (my $var = $var_fh->getline) {
        if ($var =~ /^\#/) {# header
            $out_fh->print($var);
        }
        else {
            my @columns  = split("\t", $var);
            my $ref_base = $columns[3];
            my $alt_base = $columns[4];
            next if $ref_base eq 'N';
            next if $alt_base eq '.';
            $out_fh->print($var);
        }
    }

    $var_fh->close;
    $out_fh->close;

    unless ($self->out_file) {
        unlink $var_file;
        rename $out_file, $var_file;
    }

    return 1;
}

1;
