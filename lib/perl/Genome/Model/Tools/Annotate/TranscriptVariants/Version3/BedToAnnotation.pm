package Genome::Model::Tools::Annotate::TranscriptVariants::Version3::BedToAnnotation;

use strict;
use warnings;
use Genome;
use IO::File;
use File::Temp;

class Genome::Model::Tools::Annotate::TranscriptVariants::Version3::BedToAnnotation{
    is => ['Genome::Model::Tools::Annotate'],
    has => [
        snv_file => {
            is => 'Path',
            is_input => 1,
            is_optional => 1,
            doc => 'SNV file in bed format',
        },
        indel_file => {
            is => 'Path',
            is_input => 1,
            is_optional => 1,
            doc => 'indel file in bed format',
        },
        output => {
            is => 'Path',
            is_input => 1,
            is_optional => 0,
            doc => 'File where the converted variants will be written',
        },
        extra_columns => {
            is => 'Text',
            is_input => 1,
            is_optional => 1,
            doc => "A comma delimited list of any extra columns that exist after the expected 5 in the input. Use this option if it is desired to preserve additional columns from the input file, which will then appear in output.Preserved columns must be contiguous and in order as they appear in the infile after the mandatory input columns. Any desired naming or number of columns can be specified so long as it does not exceed the actual number of columns in the file.",
        },
    ],
};

sub execute{
    my $self = shift;
    
    unless(defined $self->snv_file or defined $self->indel_file){
        $self->error_message("snv-file and/or indel-file must be defined") and die; 
    }

    my ($snv_file, $indel_file);
    $snv_file   = $self->_convert_input_file($self->snv_file) if defined $self->snv_file;
    $indel_file = $self->_convert_input_file($self->indel_file) if defined $self->indel_file;

    my $cmd = "sort -k 1,1 -k 2,2n -k 3,3n -o " . $self->output;
    if(defined $snv_file) {
        $cmd .= " " . $snv_file->filename;
    }
    if(defined $indel_file) {
        $cmd .= " " . $indel_file->filename;
    }
    if(! Genome::Sys->shellcmd(cmd => $cmd)) {
        Carp::confess "Failed to sort inputs.";
        die;
    }
    return 1;
}

sub _convert_input_file{
    my ($self, $input_file) = @_;

    my $output = File::Temp->new;
    unless ($output){
        $self->error_message("Could not open temp file, exiting") and die; 
    }

    my @columns = (($self->_variant_attributes), $self->get_extra_columns);
    my $svr = Genome::Utility::IO::SeparatedValueReader->create(
        input => $input_file,
        headers => \@columns,
        separator => "\t|\/",
        is_regex => 1,
        ignore_extra_columns => 1,
    );
    unless ($svr){
       $self->error_message("No separated value reader for $input_file, exiting") and die; 
    }
    while (my $line = $svr->next){
        #replace the '*'s with '-' as the annotator expects
        $line->{'reference'} =~ s/0|\*/-/g; 
        $line->{'variant'}   =~ s/0|\*/-/g;

        my $final;
        if($line->{'reference'} eq '-'){
            #insertions should have start == stop in Chris's format
            $final = join("\t", $line->{'chromosome'}, $line->{'start'}, $line->{'stop'} + 1, $line->{'reference'}, $line->{'variant'});
        }else{
            $final = join("\t", $line->{'chromosome'}, $line->{'start'} + 1, $line->{'stop'}, $line->{'reference'}, $line->{'variant'});
        }

        if($self->extra_columns){
            $final = join("\t", $final, map($line->{$_}, $self->get_extra_columns));
        }

        print $output $final . "\n";
    }
    #flush the buffer
    my $ofh = select $output;
    $| = 1;                         # Make $output socket hot
    print $output "";               # print nothing
    $| = 0;                         # $output socket is no longer hot
    select $ofh;
    return $output;
}

sub compare_variants{
    my ($first, $second) = @_;

    #if either parameter is undefined, return the other
    return -1 unless defined $second;
    return 1 unless defined $first;

    my ($first_chrom, $first_start, $first_stop) = split("\t", $first);
    my ($second_chrom, $second_start, $second_stop) = split("\t", $second);
    my $return_value;
    if((my $x = chr_cmp($first_chrom, $second_chrom)) != 0){
        $return_value = $x;
    }
    elsif($first_start < $second_start){
        $return_value =  -1;
    }
    elsif($second_start < $first_start){
        $return_value =  1;
    }
    elsif($first_stop < $second_stop){
        $return_value = -1;
    }
    elsif($second_stop < $first_stop){
        $return_value = 1;
    }
    else{
        $return_value = 0;
    }

    return $return_value;
}

sub chr_cmp {
    my ($a,$b) = @_;
    no warnings;
    if($a > 0 && $b > 0){
        return ($a <=> $b);
    }
    elsif($a lt $b){
        return -1;
    }
    elsif($b lt $a){
        return 1;
    }
    else{
        return 0;
    }
}

sub _variant_attributes {
    return (qw/ chromosome start stop reference variant /);
}

sub get_extra_columns {
    my $self = shift;

    my $unparsed_columns = $self->extra_columns;
    return unless $unparsed_columns;

    my @columns = split(",", $unparsed_columns);
    chomp @columns;

    return @columns;
}

1;
