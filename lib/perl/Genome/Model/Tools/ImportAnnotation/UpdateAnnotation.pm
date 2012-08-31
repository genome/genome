package Genome::Model::Tools::ImportAnnotation::UpdateAnnotation;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Model::Tools::ImportAnnotation::UpdateAnnotation {
    is  => 'Command',
    has => [
        ensembl_db => {
                       is  => "Text",
                       doc => "",
                       is_optional => 1,
                      },
        genbank_db => {
                       is  => "Text",
                       doc => "",
                       is_optional => 1,
                      },
        version    => {
                       is  => "Text",
                       doc => "annotation version string (eg 52_36n)",
                      },
        build_control => {
                           is      => "Boolean", # or something else
                           doc     => "",
                           default => 0,
                          },
    ],

};

sub sub_command_sort_position {12}

sub help_brief
{
    "Tool for getting the ensembl, genbank and combined annotation imported annotations";
}

sub help_synopsis
{
return <<EOS
fill me out
EOS
}

sub help_detail
{
    return <<EOS
This command is used for creating the ensembl, genbank and combined imported
annotation data source files, and setting up the files properly.
EOS
}

sub execute
{
    my $self = shift;

    $self->status_message("This tool is still being developed.\n come back later");
    return 1;

    my $flatfile = $self->derive_flatfile_location();
    my $gbff = "/gscmnt/sata835/info/medseq/annotation_data/human.rna.gbff"; # hard code this for now

    my $ensembl_host = "mysql1"; # or mysql2?
    my $ensembl_user = "mse";

    my $outputdir_base = "/gscmnt/sata835/info/medseq/annotation_data/";

    my $vers_number = $self->version;
    

    my $ensdir = $outputdir_base."annotation_data.ensembl_v".$vers_number;
    my $gbdir  = $outputdir_base."annotation_data.genbank_v".$vers_number;
    # they look like:
    # /gscmnt/sata835/info/medseq/annotation_data/annotation_data.something_v1/

=head1 NOTES

below from the wiki, 
B<https://gscweb.gsc.wustl.edu/wiki/User:Adukes/Imported_Annotation_Genome_Models>
 ens_dir_loc = derived location to store ensembl annotation dir
 G:M:T:ImportedAnnotation:Ensembl->create(ensembl-db=>", 
                                      output-director=>ens_dir_loc),execute
 G:M:C:B->create(model_id=><ensembl_model_id>, 
                  version=><version>, 
                  annotation_directory=>ens_dir_loc),execute 
 gen_dir_loc = derived location to store genbank annotation dir
 G:M:T:ImportedAnnotation:Genbank->create(genbank-db=>", 
                          output-director=>ens_dir_loc),execute
 G:M:C:B->create(model_id=><genbank_model_id>, 
                  version=><version>, 
                  annotation_directory=>gen_dir_loc) 
 G:M:C:B->create(model_id=><combined-annotation_model_id>, 
                  version=><version>, 
       possible_4th_build_arg_telling_me_i'm_a_composite_model?=>1),execute 

=cut

    return 1;
}

sub derive_flatfile_location
{
    my $self = shift;
    my $version = $self->version;

    my $flatfile = undef;
    $flatfile = "/gsc/var/lib/import/entrez/".$version."/Homo_sapiens.agc";
    unless(-f $flatfile)
    {
        croak "ERROR: $flatfile doesn't exist!\n";
    }
    
    return $flatfile;
}


1;
