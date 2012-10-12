package Genome::Model::Tools::Germline::Qqplot;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use Carp qw/confess/;
use POSIX qw( WIFEXITED );

class Genome::Model::Tools::Germline::Qqplot {
    is => 'Command::V2',
    has_input => [
        input_file => {
            is => 'FilesystemPath',
            doc => 'file to generate qqplot from',
            shell_args_position => 1,
        },
        output_file_basename => {
            is => 'Text',
            doc => 'basename of output file. The appropriate extension will be appended based on the image type',
            is_output => 1,
        },
    ],
    has_optional_input => [
        pvalue_column => {
            is => 'Text',
            doc => 'The column name or index of the column containing the p-values',
            default => 'p.value',
        },
        separator => {
            is => 'Text',
            doc => 'The separator between fields in the input file',
            default => '\t',
        },
        header => {
            is => 'Boolean',
            doc => 'Whether or not the input file has a header',
            default => 1,
        },
        title => {
            is => 'Text',
            doc => 'Optional title of the graph',
            default => 'Quantile-quantile plot of p-values',
        },
        image_type => {
            is => 'Text',
            doc => 'type of image to graph',
            default => 'pdf',
            valid_values => [ 'pdf', 'png' ],
        },

    ],
};

sub help_brief {
    "Graph a qqplot from a table"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt germline qqplot --input-file clinical_correlation.txt --output-file clinical_correlation_qqplot
EOS
}

sub execute {
    my $self = shift;
    $DB::single = 1;

    my $input_file = $self->input_file;

    unless(-s $input_file) {
        die $self->error_message("$input_file does not exist");
    }

        my $rscript = $self->rscript;
        my $cmd = qq{echo '$rscript'};
        $cmd .= "| R --slave";
        print "R:\n$cmd\n";
        WIFEXITED(system $cmd) or confess "Couldn't run: $cmd ($?)";

    return 1;
}

sub rscript {
    my ($self) = @_;
    my $device_setup = $self->_set_up_r_device;
    my $input_file = $self->input_file;
    my $sep = $self->separator;
    my $header = $self->header ? 'T' : 'F';
    my $title = $self->title;

    my $p_value_expression = $self->_generate_pvalue_expression;

    return <<"__RSCRIPT__";
#This code derived from http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
library(ggplot2)
$device_setup
read.table("$input_file",sep="$sep",header=$header)->x
p = x[$p_value_expression]
p = sort(p[which(!is.na(p))]) #remove any NAs
p_quantiles = sort(qchisq(p, 1, low = FALSE))
expected_p = sort(ppoints(p))
expected_p_quantiles = sort(qchisq(1 - expected_p, 1))
s <- summary(lm(p_quantiles ~ 0 + expected_p_quantiles))\$coeff
o = -log10( p )
e = -log10( expected_p )
plot = qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o)), geom="point",colour=I("blue")) + stat_abline(intercept=0,slope=1)
plot = plot + opts(title="$title")
plot = plot + scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
plot = plot + scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
plot = plot + annotate("text",label=paste(sep="","list(lambda==",signif(s[1,1],5),",s.e.==",signif(s[1, 2],5),")"), x = max(e), y = 0, hjust=1, vjust=1, parse=T)
plot
dev.off()
__RSCRIPT__

}

sub _set_up_r_device {
    my ($self) = @_;

    my $type = $self->image_type;
    my $output_file = $self->output_file_basename;
    my $dev_code;
    if($type eq 'pdf') {
        $dev_code = qq{pdf(file="$output_file.$type", height=6.5, width=9)};
    }
    elsif($type eq 'png') { 
        $dev_code = qq{png(file="$output_file.$type", height = 6.5, width = 9, units="in", res=300)};
    }
    else { 
        die "Unhandled image type " . $self->image_type;
    }
    return $dev_code;
}

sub _generate_pvalue_expression {
    my ($self) = @_;

    my $pvalue_column = $self->pvalue_column;

    if($pvalue_column =~ /^\d+$/) {
        #if it's all digits then we'll assume that it's a column index
        unless($pvalue_column > 0) {
            die $self->error_message("Column indices in R start at 1");
        }
        #not handling out of range requests. The user is responsible for handling that
        return ",$pvalue_column";
    }
    else {
        #assume it's the name of the column
        return qq{,make.names("$pvalue_column")};
    }
}

1;
