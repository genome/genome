package Genome::DruggableGene::Command::Import::Base;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::Command::Import::Base {
    is => 'Command::V2',
    is_abstract => 1,
    has => [
        version => {
            is => 'Text',
            is_input => 1,
            doc => 'Version identifier for the infile (ex 3)',
        },
        citation_base_url => {
            is => 'Text',
            is_input => 1,
            doc => 'base url string for the citation object',
        },
        citation_site_url => {
            is => 'Text',
            is_input => 1,
            doc => 'site url for the citation object',
        },
        citation_text => {
            is => 'Text',
            is_input => 1,
            doc => 'citation text for this datasource',
        },
        tmp_dir => {
            is => 'Path',
            default => '/tmp/',
            doc => 'Directory where temp files will be created',
        },        
    ],
    doc => 'Base class for importing datasets into DGI:DB',
};

sub _get_uniprot_entrez_mapping {
    my $self = shift;
    #Get mapping of Uniprot Accessions to Entrez IDs, etc
    #These can be obtained from here:
    #ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ 
    #ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz
    print "\nAttempting download of UniProt mapping file\n";
    my $mapping_file_url="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz";
    my $mapping_file_name="HUMAN_9606_idmapping_selected.tab.gz";
    my $mapping_file_path = $self->_download_file('-url'=>$mapping_file_url, '-file_name'=>$mapping_file_name);

    print "\nParsing Uniprot mapping file\n";
    my %UniProtMapping;
    open (MAPPING, $mapping_file_path) or die "can't open $mapping_file_path\n";
    while (<MAPPING>){
      my @data=split("\t",$_);
      my $uniprot_acc=$data[0];
      my $uniprot_id=$data[1];
      my $entrez_id=$data[2];
      my $ensembl_id=$data[19];
      unless ($uniprot_id){$uniprot_id="N/A";}
      unless ($entrez_id){$entrez_id="N/A";}
      unless ($ensembl_id){$ensembl_id="N/A";}
      $UniProtMapping{$uniprot_acc}{uniprot_acc}=$uniprot_acc;
      $UniProtMapping{$uniprot_acc}{entrez_id}=$entrez_id;
      $UniProtMapping{$uniprot_acc}{ensembl_id}=$ensembl_id;
    }
    close MAPPING;
    return(\%UniProtMapping);
}

sub _download_file {
    my $self = shift;
    my %args = @_;
    my $url = $args{'-url'};
    my $targetfilename;
    if ($args{'-file_name'}){
      $targetfilename = $args{'-file_name'};
    }elsif ($url=~/http.+\/(\S+)$/){ #Grab non-whitespace content after last slash to use for temp file name
      $targetfilename=$1;
    }else{
      die "could not determine file name from $url";
    }
    my $tempdir = $self->tmp_dir;
    my $targetfilepath="$tempdir"."$targetfilename";
    my $wget_cmd = "wget $url -O $targetfilepath";
    my $retval = Genome::Sys->shellcmd(cmd=>$wget_cmd);
    unless ($retval == 1){
      self->error_message('Failed to wget the specified URL');
      return;
    }
    #unzip if necessary
    if ($targetfilepath=~/\.gz$/){
      my $gunzip_cmd = "gunzip -f $targetfilepath";
      my $retval2 = Genome::Sys->shellcmd(cmd=>$gunzip_cmd);
      unless ($retval2 == 1){
        self->error_message('Failed to gunzip the specified file');
        return;
      }
      $targetfilepath=~s/\.gz$//;
    }
    print "Downloaded $targetfilepath\n";
    return $targetfilepath;
}


1;
