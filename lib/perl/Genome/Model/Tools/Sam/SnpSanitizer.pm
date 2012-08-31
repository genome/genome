package Genome::Model::Tools::Sam::SnpSanitizer;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;


class Genome::Model::Tools::Sam::SnpSanitizer {
    is  => 'Command',
    has => [
    snp_file => {
        is  => 'String',
        doc => 'The input sam/bam snp file',
    },
    ],
    has_optional => [
    out_file => {
        is  => 'String',
        doc => 'snp output file after sanitizing, default is using the same snp_file name',
    },
    ],
};


sub help_brief {
    'Sanitize samtools-pileup snp output';
}

sub help_detail {
    return <<EOS
    Sanitize samtools-pileup snp output. For samtools version r320wu1, pileup "-v" option will 
    output both snp and indel lines, and indel lines are double print(see below example). We 
    need sanitize the output to be a clean snp output to be used for all snp-related analysis. 
    This really sucks. Version r301wu1 has option "-S" to get clean SNP output. Hope "-v" option 
    can be fixed later in distribution, which would make this module obsolete.

    1       991103  C       Y       3       3       54      6       tt,,,,  ?B\@A>5
    1       991173  T       A       0       0       60      2       a\$,     BB
    1       1018280 C       C       30      0       60      1       ,-1a    ?
    1       1018280 *       -A/-A   40      0       60      1       -A      *   1   0   0
    1       1019668 A       G       16      16      60      1       G       1
    1       1019673 A       G       7       7       60      1       G       (
    1       1019689 G       T       4       4       60      1       T       # 

EOS
}


sub execute {
    my $self = shift;
    my $snp_file = $self->snp_file;

    unless (-s $snp_file) {
        $self->error_message('Can not find valid SAM snp file: '.$snp_file);
        return;
    }

    my $out_file;

    if ($self->out_file) {
        $out_file = $self->out_file;
    }
    else {
        my ($name, $path) = fileparse $snp_file;
        $out_file = $path.$name.'.sanitize';
    }


    my $out_fh = Genome::Sys->open_file_for_writing($out_file) or return;
    my $snp_fh = Genome::Sys->open_file_for_reading($snp_file) or return;

    while (my $snp = $snp_fh->getline) {
        my @columns = split("\t", $snp);
        my $ref_base = $columns[2];
        my $con_base = $columns[3];
        next if $ref_base eq $con_base;
        next if $ref_base eq '*' or $con_base =~ /\//;
        $out_fh->print($snp);
    }

    $snp_fh->close;
    $out_fh->close;

    unless ($self->out_file) {
        unlink $snp_file;
        rename $out_file, $snp_file;
    }



    return 1;
}


1;
