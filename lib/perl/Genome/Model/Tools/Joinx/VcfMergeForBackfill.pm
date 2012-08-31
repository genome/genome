package Genome::Model::Tools::Joinx::VcfMergeForBackfill;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::Tools::Joinx::VcfMergeForBackfill {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_files => {
            is => 'Text',
            is_many => 1,
            doc => 'List of bed files to sort',
            shell_args_position => 1,
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
        joinx_bin_path => {
            is => 'Text',
            doc => 'path to the joinx binary to use. This tool is being released before joinx vcf-merge will be released. This will go away when it is.',
        },
    ],
};

sub help_brief {
    "Sorts one or more bed files."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt joinx vcf-merge a.vcf b.vcf ... --output-file merged.vcf
EOS
}

sub execute {
    my $self = shift;
    if(defined($self->use_bgzip) && not defined($self->output_file)){
       die $self->error_message("If use_bgzip is set, output_file must also be set, otherwise binary nonsense will spew forth."); 
    }
    my $output = "-";
    $output = $self->output_file if (defined $self->output_file);

    # Grep out empty files
    my @inputs = grep { -s $_ } $self->input_files;

    # If all input files are empty, make sure the output file at least exists
    unless (@inputs) {
        if (defined $self->output_file) {
            unless (system("touch $output") == 0) {
                die $self->error_message("Failed to touch $output");
            }
        }
        return 1;
    }
    
    if($self->use_bgzip){
        my @new_inputs;
        for my $input (@inputs){
            push @new_inputs, "<(zcat $input)";
        }
        @inputs = @new_inputs;
    }

    unless($self->joinx_bin_path){
        $self->joinx_bin_path("joinx");
    }

    #compose the command, using awk to compose bed compatible records
    my $cmd = $self->joinx_bin_path . " vcf-merge " . join(" ", @inputs);
    if(defined($self->output_file) && not defined($self->use_bgzip)){
        $cmd .= " -o $output" if defined($self->output_file);
    } elsif ( defined($self->use_bgzip) && defined($self->output_file) ){
        $cmd .= " | grep -v ^# | cut -f1-5 | awk \'BEGIN { OFS=\\\"\\t\\\"; }  { print \\\$1,\\\$2,\\\$2,\\\$5; }\' | bgzip -c > $output";
        $cmd = "bash -c \"$cmd\"";
    }
        
    my %params = (
        cmd => $cmd,
        allow_zero_size_output_files=>1,
    );
    $params{output_files} = [$output] if $output ne "-";
    Genome::Sys->shellcmd(%params);

    return 1;
}

1;
