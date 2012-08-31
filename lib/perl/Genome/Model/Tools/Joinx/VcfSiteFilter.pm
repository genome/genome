package Genome::Model::Tools::Joinx::VcfSiteFilter;

use strict;
use warnings;

use Genome;
use Data::Dumper;

our $MINIMUM_JOINX_VERSION = 1.3; #Not present in earlier versions

class Genome::Model::Tools::Joinx::VcfSiteFilter {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Vcf File to filter',
            shell_args_position => 1,
        },
        min_fail_filter => {
            is => 'Number',
            doc => 'Minimum fraction of failed genotypes to fail a site',
        },
    ],
    has_optional_input => [
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => 'The output file (defaults to stdout)',
        },
        use_bgzip => {
            is => 'Boolean',
            doc => 'zcats the input files into stdin, and bgzips the output',
            default => 0,
        },
    ],
};

sub help_brief {
    "Filters out sites in a VCF file based on the number of samples that had genotypes that were filtered out."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt joinx vcf-site-filter --min-fail-filter 0.5 a.vcf ... --output-file merged.vcf
EOS
}

sub execute {
    my $self = shift;

    if($self->use_version < $MINIMUM_JOINX_VERSION) {
        die $self->error_message("This module requires joinx version $MINIMUM_JOINX_VERSION or higher to function correctly.");
    }
    
    if(defined($self->use_bgzip) && not defined($self->output_file)){
       die $self->error_message("If use_bgzip is set, output_file must also be set, otherwise binary nonsense will spew forth."); 
    }
    my $output = "-";
    $output = $self->output_file if (defined $self->output_file);

    my $input_file = $self->input_file;

    unless(-s $input_file) {
        die $self->error_message("Unable to read from $input_file");
    }

    if($self->use_bgzip){
        $input_file = "<(zcat $input_file)";
    }

    #compose the command, using awk to compose bed compatible records
    my $cmd = $self->joinx_path . " vcf-site-filter --min-fail-filter " . $self->min_fail_filter . " --input-file $input_file";
    if(defined($self->output_file) && not defined($self->use_bgzip)){
        $cmd .= " -o $output" if defined($self->output_file);
    } elsif ( defined($self->use_bgzip) && defined($self->output_file) ){
        $cmd .= " | bgzip -c > $output";
        $cmd = "bash -c \"$cmd\"";
    }
    $self->status_message($cmd);
        
    my %params = (
        cmd => $cmd,
    );
    $params{output_files} = [$output] if $output ne "-";
    $params{skip_if_output_is_present} = 0;
    Genome::Sys->shellcmd(%params);

    return 1;
}

1;
