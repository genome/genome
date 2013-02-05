package Genome::Model::Tools::Germline::Manhattan;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use Carp qw/confess/;
use POSIX qw( WIFEXITED );

class Genome::Model::Tools::Germline::Manhattan {
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
        image_type => {
            is => 'Text',
            doc => 'type of image to graph',
            default => 'pdf',
            valid_values => [ 'pdf', 'png' ],
        },
        input_is_clinical_correlation => {
            is => 'Boolean',
            doc => 'invoke special processing for clinical correlation output',
            default => 1,
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

    my $handle_file = $self->input_is_clinical_correlation ? $self->_clinical_correlation_handling("data") : $self->_other_file_handling("data");
    my $drawing_lib = __FILE__ . ".R";

    return <<"__RSCRIPT__";
#This code derived from http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
library(ggplot2)
$device_setup
source("$drawing_lib")
$handle_file
manhattan(data)
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

sub _clinical_correlation_handling {
    my ($self, $variable_name) = @_;
    #here we are assuming that all is known. Probably suboptimal and more than a little redundant with generate_pvalue_expression etc
    my $input = $self->input_file;
    my $sep = $self->separator;
    my $pvalue_expression = $self->_generate_pvalue_expression;

    my $clinical_correlation_handling = <<"__RCODE__";
$variable_name<-read.table("$input",header=T,sep="$sep")
strs=strsplit(as.character($variable_name\$x),"_")
$variable_name=data.frame(CHR=sapply(strs,"[[",1),BP=as.numeric(sapply(strs,"[[",2)), P=${variable_name}[$pvalue_expression])
$variable_name\$CHR = sub("^X(.+)",replacement="\\\\\\\\1",$variable_name\$CHR,fixed=FALSE)
$variable_name\$CHR = sub("^X",replacement="23",$variable_name\$CHR,fixed=FALSE)
$variable_name\$CHR = as.numeric(sub("^Y",replacement="24",$variable_name\$CHR,fixed=FALSE))
__RCODE__
    
    return $clinical_correlation_handling;
}

sub _other_file_handling {
    my ($self, $variable_name) = @_;

    my $input = $self->input_file;
    my $sep = $self->separator;
    my $pvalue_expression = $self->_generate_pvalue_expression;
    my $header = $self->header ? 'T' : 'F';

    #NOTE this assumes everything is generally peachy keen with the exception
    #of the p-value expression stuff.  I assume this will likely be used only
    #for raw PLINK output but it may break for other inputs
    my $other_file_handling = <<"__RCODE__";
$variable_name<-read.table("$input",header=$header,sep="$sep")
$variable_name=transform($variable_name,P=${variable_name}[$pvalue_expression])
__RCODE__

    return $other_file_handling;
}


1;
