package Genome::Model::Tools::Picard::PlotRnaSeqMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::PlotRnaSeqMetrics {
    is  => 'Command::V2',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The metrics file output from CollectRnaSeqMetrics.',
        },
        output_file => {
            is => 'String',
            doc => 'The output plot.',
        },
        label => {
            is_optional => 1,
            is => 'String',
            doc => 'A label to add to the plot title.',
        },
    ],
};

sub help_brief {
    'Program to plot metrics about the alignment of RNA to various functional classes of loci in the genome: coding, intronic, UTR, intergenic, ribosomal.';
}

sub help_detail {
    return <<EOS
    The command is based on the output of CollectRnaSeqMetrics,
    for Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CollectRnaSeqMetrics
EOS
}

sub execute {
    my $self = shift;

    my $metrics_hashref = Genome::Model::Tools::Picard::CollectRnaSeqMetrics->parse_file_into_metrics_hashref($self->input_file);
    #print Data::Dumper::Dumper($metrics_hashref);
    my $utr_bases = $metrics_hashref->{UTR_BASES};
    my $coding_bases = $metrics_hashref->{CODING_BASES};
    my $ribosomal_bases = $metrics_hashref->{RIBOSOMAL_BASES};
    my $intergenic_bases = $metrics_hashref->{INTERGENIC_BASES};
    my $intronic_bases = $metrics_hashref->{INTRONIC_BASES};

    my $title = 'RNASeq Alignment Breakdown';
    if ($self->label) {
        $title .= ' : '. $self->label;
    }
    my ($tmp_fh,$tmp_file) = Genome::Sys->create_temp_file;
    print $tmp_fh 'library(plotrix)' ."\n";
    print $tmp_fh 'slices <- c('.$utr_bases .', '. $coding_bases .', '. $ribosomal_bases .', '. $intronic_bases .', '. $intergenic_bases .')' ."\n";
    print $tmp_fh 'lbls <- c("UTR", "CODING", "RIBOSOMAL", "INTRONIC", "INTERGENIC")' ."\n";
    print $tmp_fh 'pct <- round(slices/sum(slices)*100)' ."\n";
    print $tmp_fh 'lbls <- paste(lbls, pct)' ."\n"; # add percents to labels
    print $tmp_fh 'lbls <- paste(lbls,"%",sep="")' ."\n"; # ad % to labels 
    print $tmp_fh 'png("'.$self->output_file .'", bg="transparent", width=800, height=600)' ."\n";
    print $tmp_fh 'pie3D(slices,labels=lbls,explode=0.1,main="'. $title .'")' ."\n";
    print $tmp_fh 'dev.off()' ."\n";
    $tmp_fh->close;
    my $cmd = 'Rscript '. $tmp_file;
    eval {
        Genome::Sys->shellcmd(
            cmd => $cmd,
        );
    };
    if ($@) {
        warn('The Rscript did not run correctly: '. $@);
    }
    return 1;
}

1;
