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

sub _create_citation {
    my $self = shift;
    my $source_db_name = shift;
    my $source_db_version = shift;
    my $base_url = shift;
    my $site_url = shift;
    my $citation_text = shift;
    my $full_name = shift;
    return Genome::DruggableGene::Citation->create(
        source_db_name => $source_db_name,
        source_db_version => $source_db_version,
        base_url => $base_url,
        site_url => $site_url,
        citation => $citation_text,
        full_name => $full_name,
    );
}

sub _create_drug_name_report {
    my $self = shift;
    my ($name, $citation, $nomenclature, $description) = @_;
    my %params = (
        name => uc $name,
        nomenclature => $nomenclature,
        citation => $citation,
        description => $description,
    );

    my $drug_name_report = Genome::DruggableGene::DrugNameReport->get(%params);
    return $drug_name_report if $drug_name_report;
    return Genome::DruggableGene::DrugNameReport->create(%params);
}

sub _create_drug_alternate_name_report {
    my $self = shift;
    my ($drug_name_report, $alternate_name, $nomenclature, $description) = @_;
    my %params = (
        drug_id => $drug_name_report->id,
        alternate_name => uc $alternate_name,
        nomenclature => $nomenclature,
        description => $description,
    );

    my $drug_alternate_name_report = Genome::DruggableGene::DrugAlternateNameReport->get(%params);
    return $drug_alternate_name_report if $drug_alternate_name_report;
    return Genome::DruggableGene::DrugAlternateNameReport->create(%params);
}

sub _create_drug_category_report {
    my $self = shift;
    my ($drug_name_report, $category_name, $category_value, $description) = @_;
    my %params = (
        drug_id => $drug_name_report->id,
        category_name => $category_name,
        category_value => lc $category_value,
        description => $description,
    );
    my $drug_category_report = Genome::DruggableGene::DrugCategoryReport->get(%params);
    return $drug_category_report if $drug_category_report;
    return Genome::DruggableGene::DrugCategoryReport->create(%params);
}

sub _create_gene_name_report {
    my $self = shift;
    my ($name, $citation, $nomenclature, $description) = @_;
    my %params = (
        name => uc $name,
        nomenclature => $nomenclature,
        citation => $citation,
        description => $description,
    );

    if($name ne 'NA'){
        my $gene_name_report = Genome::DruggableGene::GeneNameReport->get(%params);
        return $gene_name_report if $gene_name_report;
    }
    return Genome::DruggableGene::GeneNameReport->create(%params);
}

sub _create_gene_alternate_name_report {
    my $self = shift;
    my ($gene_name_report, $alternate_name, $nomenclature, $description, $case) = @_;
    if ($case eq 'upper'){$alternate_name = uc $alternate_name}
    my %params = (
        gene_id => $gene_name_report->id,
        alternate_name => $alternate_name,
        nomenclature => $nomenclature,
        description => $description,
    );
    my $gene_alternate_name_report = Genome::DruggableGene::GeneAlternateNameReport->get(%params);
    return $gene_alternate_name_report if $gene_alternate_name_report;
    return Genome::DruggableGene::GeneAlternateNameReport->create(%params);
}

sub _create_gene_category_report {
    my $self = shift;
    my ($gene_name_report, $category_name, $category_value, $description) = @_;
    my %params = (
        gene_id => $gene_name_report->id,
        category_name => $category_name,
        category_value => $category_value,
        description => $description,
    );
    my $gene_category_report = Genome::DruggableGene::GeneCategoryReport->get(%params);
    return $gene_category_report if $gene_category_report;
    return Genome::DruggableGene::GeneCategoryReport->create(%params);
}

sub _create_interaction_report {
    my $self = shift;
    my ($citation, $drug_name_report, $gene_name_report, $description) = @_;
    my %params = (
        gene_id => $gene_name_report->id,
        drug_id => $drug_name_report->id,
        citation => $citation,
        description =>  $description,
    );

    my $interaction = Genome::DruggableGene::DrugGeneInteractionReport->get(%params);
    return $interaction if $interaction;
    return Genome::DruggableGene::DrugGeneInteractionReport->create(%params);
}

sub _create_interaction_report_attribute {
    my $self = shift;
    my ($interaction, $name, $value) = @_;
    my %params = (
        drug_gene_interaction_report => $interaction,
        name => $name,
        value => lc $value,
    );
    my $attribute = Genome::DruggableGene::DrugGeneInteractionReportAttribute->get(%params);
    return $attribute if $attribute;
    return Genome::DruggableGene::DrugGeneInteractionReportAttribute->create(%params);
}

sub _destroy_and_rebuild_pubchem_and_drug_groups {
    my $citation_text = Genome::DruggableGene::Citation->get(source_db_name => 'PubChem')->citation;
    #nuke drug groups
    Genome::DruggableGene::Command::DrugNameGroup::RemoveAll->execute();
    #nuke pubchem
    Genome::DruggableGene::Command::RemoveDatasource->execute(source_db_name => 'PubChem');
    #reimport pubchem
    Genome::DruggableGene::Command::Import::Pubchem->execute(postgres_host => 'postgres', scratch_dir => '/tmp/', citation_text => $citation_text, generate_uuids_locally => 1); #TODO: handle different postgres hosts
    #regroup_drugs
    Genome::DruggableGene::Command::DrugNameGroup::Generate->execute();
    return 1;
}

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
