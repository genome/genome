package Genome::Model::Tools::Far::Trimmer;

use strict;
use warnings;

use Genome;

my $VALID_FORMATS = ['fastq','fastq-illumina13','fastq-sanger','fastq-solexa','fasta','csfastq','csfasta'];

class Genome::Model::Tools::Far::Trimmer {
    is => ['Genome::Model::Tools::Far::Base'],
    has_input => [
        source => {
            is => 'Text',
            doc => 'Input file containing reads to be trimmed',
        },
        target => {
            is => 'Text',
            doc => 'Output file containing trimmed reads',
        },
    ],
    has_optional_input => [
        source2 => {
            is => 'Text',
            doc => 'Second input file containing reads to be trimmed',
        },
        params => {
            is => 'Text',
            doc => 'A list of params to be APPENDED to the command line.',
        },
        adapter => {
            is => 'Text',
            doc => 'String (adaptor sequence) to be removed',
        },
        format => {
            is => 'Text',
            doc => 'input file format, output will be in the same format',
            valid_values => $VALID_FORMATS,
        },
        cut_off => {
            is => 'Text',
            doc => 'Max. nr of allowed mismatches + indels for alignment per 10 bases (far default = 2)',
        },
        min_readlength => {
            is => 'Text',
            doc => 'minimum readlength in basepairs after adapter removal - read will be discarded otherwise. far default:18',
        },
        max_uncalled => {
            is => 'Text',
            doc => 'nr of allowed uncalled bases in a read. far default: 0',
        },
        min_overlap => {
            is => 'Text',
            doc => 'minimum required overlap of adapter and sequence in basepairs. far default: 10',
        },
        nr_threads => {
            is            => 'Text',
            doc           => 'Number of threads to use. far default: 1',
        },
        trim_end => {
            is          => 'Text',
            doc         => 'Decides on which end adapter removal is performed. far default: right',
            valid_values => ['right','left','any','left_tail','right_tail'],
        },
        far_output => {
            is => 'Text',
            doc => 'Redirect the stdout from the far command.',
        },
        trim_reverse_complement => {
            is => 'Boolean',
            doc => 'Setting this flag will also trim the reverse complement of the --adapter argument',
            default => 0,
        },
    ],
    has => [
        _far_cmd => {
            is => 'Text', 
            doc => 'Far command to run',
            is_optional => 1,
        },
        threads => {
            is            => 'Text',
            doc           => 'Deprecated, use --nr-threads',
            is_optional => 1,
        },
        min_read_length => {
            is => 'Text',
            doc => 'Deprecated, use --min-readlength',
            is_optional => 1,
        },
        file_format => {
            is => 'Text',
            doc => 'Deprecated, use --format',
            valid_values => $VALID_FORMATS,
            is_optional => 1,
        },
        adaptor_sequence => {
            is => 'Text',
            doc => 'Deprecated, use --adapter',
            is_optional => 1,
        },
    ],
};


sub execute {
    my $self = shift;
    my $far_cmd = $self->far_path;

    #extract adapter and format from params string if supplied and set it to the adapter property.  this allows us to create an adapters.fa file if trim_reverse_complement is set and set format to default fastq
    my @cmd_line_params = split(/\s+/, $self->params);
    my @new_params;
    for (my $i=0; $i < scalar @cmd_line_params; $i++){
        my $param = $cmd_line_params[$i];   
        if ($param eq '-as' or $param eq '--adapter'){
            $i++;
            my $sequence = $cmd_line_params[$i];
            $self->adapter($sequence);
        }elsif ($param eq '-f' or $param eq '--format'){
            $i++;
            my $format = $cmd_line_params[$i];
            if (defined $self->format){
                die $self->error_message("format specified in params string and also provided as property, only use one!");
            }
            $self->format($format);
        }else{
            push @new_params, $param;
        }
    }
    $self->params(join(" ", @new_params));

    #make sure adapter is uppercase and that a format is set
    $self->adapter(uc($self->adapter));
    $self->format("fastq") unless defined $self->format;


    #check deprecated param names and use the updated names(which are the same as the corresponding far param)  Eventually we will die here, and then remove the options
    for (['threads', 'nr_threads'], ['min_read_length', 'min_readlength'], ['adaptor_sequence', 'adapter'], ['file_format', 'format']){
        my ($deprecated, $new) = @$_;
        if ( defined ($self->$deprecated) ){
            $self->warning_message("Property $deprecated is deprecated!  Please use property $new!");
            if ( defined ($self->$new) ){
                die $self->error_message("Deprecated property $deprecated and corresponding propery $new are both defined, define only property $new!");
            }
            $self->$new($self->$deprecated);
        }
    }

    #create far command, check for duplicate args in obj properties and the params string
    my @params = (qw/source target source2 format nr_threads trim_end min_readlength max_uncalled min_overlap/);
    for my $param (@params){
        next unless defined $self->$param;
        my $cmd_param = $param;
        $cmd_param =~ s/_/-/g;
        $cmd_param = '--'.$cmd_param;
        if (defined ($self->params) and $self->params =~ /$cmd_param/){
            die $self->error_message("Property $param is defined when it already exists in the supplied far params: ".$self->params);
        }
        $far_cmd.= " $cmd_param ".$self->$param;
    }

    #if the trim_reverse_complement flag is set, we need to 
    if ($self->trim_reverse_complement){
        my ($tmp_fh,$tmp_path) = Genome::Sys->create_temp_file;
        my $revcomp_adapter = reverse $self->adapter;
        $revcomp_adapter =~ tr/TCGA/AGCT/;
        $tmp_fh->print(">adapter1(".$self->adapter.")\n");
        $tmp_fh->print($self->adapter."\n");
        $tmp_fh->print(">adapter2($revcomp_adapter)\n");
        $tmp_fh->print("$revcomp_adapter\n");
        $tmp_fh->close;
        $far_cmd.= ' --adapters '.$tmp_path;

    }else{
        $far_cmd.= ' --adapter '.$self->adapter;
    }

    if (defined($self->params)) {
        $far_cmd .= ' '. $self->params;
    }

    if (defined($self->far_output)) {
        $far_cmd .= ' > '. $self->far_output;
    }

    Genome::Sys->shellcmd(
        cmd => $far_cmd,
        input_files => [$self->source],
        skip_if_output_is_present => 0,
    );
    return 1;
}



1;
