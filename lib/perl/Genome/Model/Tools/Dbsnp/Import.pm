package Genome::Model::Tools::Dbsnp::Import;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dbsnp::Import {
    is => 'Genome::Model::Tools::Dbsnp',
    has => [
        flat_file_url => {
            is => 'Text',
            is_input => 1,
            doc => 'Path to dbsnp flat files on ftp site',
        },
        filename_pattern => {
            is => 'Text',
            is_input => 1,
            default => 'ds_flat_chX.flat.gz',
            doc => 'String representing the naming scheme of the flatfiles in the input-directory.  The first instance of the character "X" will be replaced with the chromosome designator to create the file name (ex: ds_flat_chX.flat)',
        },
        output_file => {
            is => 'Path',
            is_output => 1,
            doc => 'Path to the final output file in .bed',
        },
        
    ],
    has_optional => [
        reference_coordinates => {
            is => 'String',
            doc => 'reference_coordinates whose coordinates will be used, regex syntax accepted for matching multiple   patch levels'
        },
        contig_name_translation_file => {
            is => 'Path',
            doc => 'File path that contains translations of contig names',
        },
        from_names_column => {
            is => 'Number',
            doc => '0-based column number containing names you want to translate from',
        },
        to_names_column => {
            is => 'Number',
            doc => '0-based column number containing names you want to translate to',
        },
        chromosome_names => {
            is => 'String',
            is_many => 1,
        },
    ],
};

sub help_brief {
    'Create unfiltered bed file from Dbsnp flat files'
}

sub help_synopsis {
    return <<EOS
gmt dbsnp import --flat_file_url ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/ASN1_flat/ --output_file output.bed 
EOS
}

sub help_detail {
    return <<EOS
This command takes a directory of Dbsnp flat files, parses them out, and creates a single, sorted bed file suitable
for creating a Genome::Model::Build::ImportedVariationList
EOS
}

sub execute {
    my $self = shift;

    my $temp_dir = Genome::Sys->base_temp_directory();

    my $output_file = join("/", $temp_dir, "/unsorted.bed");
    $self->debug_message("Using reference coordinates ".$self->reference_coordinates);

    for my $chromosome ($self->chromosome_names){
        my $flatfile = $self->filename_pattern;
        $flatfile =~ s/X/$chromosome/;
        my $file_url = join('/', $self->flat_file_url, $flatfile);
        my %params = (flatfile => $file_url, output_file => $output_file);
        if ($self->reference_coordinates) {
            $params{reference_coordinates} = $self->reference_coordinates;
        }
        my $cmd = Genome::Model::Tools::Dbsnp::Import::Flatfile->create(%params);
        unless($cmd->execute){
            $self->error_message("Failed to import flatfile $flatfile: $@");
            return 0;
        }
    }
    my $flatfile = $self->filename_pattern;
    $flatfile =~ s/X/Un/;
    my $file_url = join('/', $self->flat_file_url, $flatfile);
    my $cmd = Genome::Model::Tools::Dbsnp::Import::Flatfile->create(
                ($self->reference_coordinates ? (reference_coordinates => $self->reference_coordinates):()),
                flatfile => $file_url, 
                output_file => $output_file,
                ($self->contig_name_translation_file ? (contig_name_translation_file => $self->contig_name_translation_file) : ()),
                ($self->from_names_column ? (from_names_column => $self->from_names_column) :()),
                ($self->to_names_column ? (to_names_column => $self->to_names_column) : ()),
                use_contig => 1);
    unless($cmd->execute) {
        $self->error_message("Failed to import flatfile $flatfile: $@");
        return 0;
    }

    my @output_files = ($output_file);

    unless(Genome::Model::Tools::Joinx::Sort->execute(input_files => \@output_files, output_file => $self->output_file)){
        $self->error_message("Failed to merge and sort imported flatfiles: $@");
        return 0;
    }
    #TODO: do gabe's white/black listing, make a feature list out of the filtered bed file and use it to create a new build of the dbsnp model

    return 1;
}

1;
