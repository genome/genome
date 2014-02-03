package Genome::DruggableGene::Command::Import::TherapeuticTargetDatabase;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::Command::Import::TherapeuticTargetDatabase {
    is => 'Genome::DruggableGene::Command::Import::Base',
    has => [
        targets_raw_url => {
            is => 'Text',
            is_input => 1,
            default => 'http://bidd.nus.edu.sg/group/cjttd/TTD_download.txt',
            doc => 'URL PATH. URL to TTD targets information in raw format',
        },
        tmp_dir => {
            is => 'Path',
            default => '/tmp/',
            doc => 'Directory where temp files will be created',
        },
        drugs_crossmatching_url => {
            is => 'Text',
            is_input => 1,
            default => 'http://bidd.nus.edu.sg/group/cjttd/TTD_crossmatching.txt',
            doc => 'URL PATH. URL to TTD crossmatching information',
        },
        drug_synonyms_url => {
            is => 'Text',
            is_input => 1,
            default => 'http://bidd.nus.edu.sg/group/cjttd/Synonyms.txt',
            doc => 'URL PATH. URL to TTD drugs synonyms file',
        },
        interactions_outfile => {
            is => 'Path',
            is_input => 1,
            example_values => ['/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/TTD_WashU_INTERACTIONS.tsv'],
            doc => 'PATH.  Path to .tsv file for drug gene interactions',
        },
        citation_base_url => {
            default => 'http://bidd.nus.edu.sg/group/cjttd/ZFTTD',
        },
        citation_site_url => {
            default => 'http://bidd.nus.edu.sg/group/ttd/ttd.asp',
        },
        citation_text => {
            default => "Update of TTD: Therapeutic Target Database. Zhu F, Han BC, ..., Zheng CJ, Chen YZ. Nucleic Acids Res. 38(suppl 1):D787-91, 2010. PMID: 19933260.",
        },
    ],
    has_optional_transient => [
        version => {
            is => 'Text',
            default => 'this parameter is ignored by the TTD importer',
            doc => 'this parameter is ingored by the TTD importer--the version is taken from the input files',
        },
    ],
    doc => 'Parse TTD download, crossmatching and synonyms files and import interactions',
};

sub help_synopsis {
    return <<HELP
genome druggable-gene import therapeutic-target-database 
HELP
}

my %UniProtMapping;

sub execute {
    my $self = shift;
    %UniProtMapping=%{$self->_get_uniprot_entrez_mapping()}; #Load UniProt to Entrez mapping information from file (For Uniprot -> Entrez mapping)
    $self->input_to_tsv();
    return 1;
}

sub input_to_tsv {
    my $self = shift;
    my $targets_raw_url = $self->targets_raw_url;
    my $drugs_crossmatching_url = $self->drugs_crossmatching_url;
    my $drug_synonyms_url = $self->drug_synonyms_url;

    #Create interactions output file
    my $interactions_outfile = $self->interactions_outfile;
    my $interactions_fh = IO::File->new($interactions_outfile, 'w');
    my $interactions_header = join("\t", 'drug_id', 'drug_name', 'drug_synonyms', 'drug_cas_number', 'drug_pubchem_cid', 'drug_pubchem_sid', 'target_id', 'target_name', 'target_synonyms', 'target_uniprot_id', 'target_entrez_id', 'target_ensembl_id', 'interaction_type');
    $interactions_fh->print($interactions_header, "\n");

    #Get the data in order
    my $targets_path = $self->_download_file('-url'=>$targets_raw_url);
    my $crossmatch_path = $self->_download_file('-url'=>$drugs_crossmatching_url);
    my $synonyms_path = $self->_download_file('-url'=>$drug_synonyms_url);

    my ($targets, $version) = $self->_parse_targets_file($targets_path);
    $self->version($version) if $version;
    my $drugs = $self->_parse_crossmatch_file($crossmatch_path);
    $self->_parse_synonyms_file($synonyms_path, $drugs);

    #Load UniProt to Entrez mapping information from file - This will be used to obtain Entrez IDs from UniProt accessions provided in Drugbank records
    my %UniProtMapping=%{$self->_get_uniprot_entrez_mapping()};

    #Write data to the file
    for my $target_id (keys %{$targets}){
            #Target Uniprot Id
            my $target_uniprot_id =  pop @{$targets->{$target_id}{'UniProt ID'}};
            $target_uniprot_id = 'N/A' unless $target_uniprot_id;

            #Retrieve Entrez/Ensembl IDs for interaction protein (if available)
            my $entrez_id = "N/A";
            my $ensembl_id = "N/A";
            if ($UniProtMapping{$target_uniprot_id}){
              $entrez_id = $UniProtMapping{$target_uniprot_id}{entrez_id};
              $ensembl_id = $UniProtMapping{$target_uniprot_id}{ensembl_id};
            }

            #Target Name
            my $target_name = pop @{$targets->{$target_id}{'Name'}};

            #Target Synonyms
            my $target_synonyms;
            $target_synonyms = join(";", @{$targets->{$target_id}{'Synonyms'}}) if $targets->{$target_id}{'Synonyms'};
            $target_synonyms = 'N/A' unless $target_synonyms;

        for my $drug_id (keys %{$targets->{$target_id}{'drugs'}}){
            #Drug Name
            my (undef, $drug_name) = split("\t", $targets->{$target_id}{'drugs'}{$drug_id});

            #CAS Number
            my $drug_cas_number = $drugs->{$drug_id}{'CAS Number'};
            if($drug_cas_number){
                $drug_cas_number =~ s/CAS //;
            }else{
                $drug_cas_number = 'N/A';
            }

            #PubChem CID
            my $drug_pubchem_cid = $drugs->{$drug_id}{'PubChem CID'};
            if($drug_pubchem_cid){
                $drug_pubchem_cid =~ s/CID //;
            }else{
                $drug_pubchem_cid = 'N/A';
            }

            #PubChem SID
            my $drug_pubchem_sid = $drugs->{$drug_id}{'PubChem SID'};
            if($drug_pubchem_sid){
                $drug_pubchem_sid =~ s/SID //;
            }else{
                $drug_pubchem_sid = 'N/A';
            }

            #Drug Synonyms
            my $drug_synonyms = $drugs->{$drug_id}{'synonyms'};
            $drug_synonyms = "N/A" unless $drug_synonyms;

            #Interaction Type
            my $interaction_type = $self->_determine_interaction_type($targets->{$target_id}, $drug_name, $drug_id);
            $interaction_type = 'N/A' unless $interaction_type;

            $interactions_fh->print(join("\t", $drug_id, $drug_name, $drug_synonyms, $drug_cas_number, $drug_pubchem_cid, $drug_pubchem_sid, $target_id, $target_name, $target_synonyms, $target_uniprot_id, $entrez_id, $ensembl_id, $interaction_type), "\n");
        }
    }

    $interactions_fh->close;
    return 1;
}

sub _parse_targets_file {
    my $self = shift;
    my $targets_path = shift;
    my $version;
    my $targets = {};
    my $fh = IO::File->new($targets_path, 'r');

    while(my $line = <$fh>){
        next unless $line;
        chomp $line;
        $line =~ s/\r//g;
        if($line =~ m/^Version/){
            $version = $line;
            $version =~ s/Version //i;
            next;
        }elsif($line =~ m/^TTD\w+/){
            my ($id, $key, $value, @extra_fields) = split("\t", $line);
            $key = 'drugs' if $key eq 'Drug(s)';
            if(!$targets->{$id}){
                $targets->{$id} = {};
            }
            if(!$targets->{$id}{$key}){
                if($key eq 'drugs'){
                    $targets->{$id}{'key'} = {};
                }else{
                    $targets->{$id}{$key} = [];
                }
            }
            if($key eq 'drugs'){
                my $drug_id = shift @extra_fields;
                $targets->{$id}{$key}{$drug_id} = join("\t", $drug_id, $value, @extra_fields);
            }else{
                push(@{$targets->{$id}{$key}}, join("\t", $value, @extra_fields));
            }
        }else{
            #skip this line
            next;
        }
    }
    $fh->close;
    return ($targets, $version);
}

sub _parse_crossmatch_file {
    my $self = shift;
    my $crossmatch_path = shift;
    my $drugs = {};
    my $fh = IO::File->new($crossmatch_path, 'r');

    while(my $line = <$fh>){
        next unless $line;
        chomp $line;
        $line =~ s/\r//g;
        if($line =~ m/^D\w+/){
            my ($id, $pointless_number, $key, $value) = split("\t", $line);
            if(!$drugs->{$id}){
                $drugs->{$id} = {};
            }
            $drugs->{$id}{$key} = $value;
        }
    }

    $fh->close;
    return $drugs;
}

sub _parse_synonyms_file {
    my $self = shift;
    my ($synonyms_path, $drugs) = @_;
    my $fh = IO::File->new($synonyms_path, 'r');

    while(my $line = <$fh>){
        next unless $line;
        chomp $line;
        $line =~ s/\r//g;
        if($line =~ m/^D\w+/){
            my ($id, $name, $synonyms) = split("\t", $line);
            $drugs->{$id}{'synonyms'} = $synonyms;
        }
    }

    $fh->close;
    return $drugs;
}

sub _determine_interaction_type{
    my $self = shift;
    my ($target, $drug_name, $drug_id) = @_;
    my $interaction_type;

    my $drug_name_and_id_string = join("\t", $drug_name, $drug_id);
    my @interaction_types = qw(Inducer Inhibitor Stimulator Agonist Blocker Antisense Antagonist Binder Activator Antibody Modulator Cofactor); #TODO: is this complete?
    for my $type (@interaction_types){
        my $drugs_by_type = $target->{$type};
        for my $drug_info (@$drugs_by_type){
            if($drug_info eq $drug_name_and_id_string){
                $interaction_type = $type;
                $interaction_type = lc $interaction_type;
            }
        }
    }

    return $interaction_type;
}


1;
