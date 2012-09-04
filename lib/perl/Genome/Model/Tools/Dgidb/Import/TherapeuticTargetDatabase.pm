package Genome::Model::Tools::Dgidb::Import::TherapeuticTargetDatabase;

use strict;
use warnings;

use Genome;

my $high = 750000;
UR::Context->object_cache_size_highwater($high);

class Genome::Model::Tools::Dgidb::Import::TherapeuticTargetDatabase {
    is => 'Genome::Model::Tools::Dgidb::Import::Base',
    has => [
        targets_raw_url => {
            is => 'Text',
            is_input => 1,
            default => 'http://bidd.nus.edu.sg/group/cjttd/TTD_download.txt',
            doc => 'URL PATH. URL to TTD targets information in raw format',
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
            default => '/tmp/TTD_WashU_INTERACTIONS.tsv',
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

sub execute {
    my $self = shift;
    $self->input_to_tsv();
    $self->import_tsv();
    $self->_destroy_and_rebuild_pubchem_and_drug_groups();
    return 1;
}

sub import_tsv {
    my $self = shift;
    my $interactions_outfile = $self->interactions_outfile;
    $self->preload_objects;
    my @interactions = $self->import_interactions($interactions_outfile);
    return 1;
}

sub import_interactions {
    my $self = shift;
    my $interactions_outfile = shift;
    my $version = $self->version;
    my @interactions;
    my @headers = qw( drug_id drug_name drug_synonyms drug_cas_number drug_pubchem_cid drug_pubchem_sid target_id target_name target_synonyms target_uniprot_id interaction_types );
    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        input => $interactions_outfile,
        headers => \@headers,
        separator => "\t",
        is_regex => 1,
    );

    my $citation = $self->_create_citation('TTD', $version, $self->citation_base_url, $self->citation_site_url, $self->citation_text);

    $parser->next; #eat the headers
    while(my $interaction = $parser->next){
        my $drug_name = $self->_import_drug($interaction, $citation);
        my $gene_name = $self->_import_gene($interaction, $citation);
        my $drug_gene_interaction = $self->_create_interaction_report($citation, $drug_name, $gene_name, '');
        push @interactions, $drug_gene_interaction;
        my @interaction_types = split('; ', $interaction->{interaction_types});
        for my $interaction_type (@interaction_types){
            my $type_attribute = $self->_create_interaction_report_attribute($drug_gene_interaction, 'interaction_type', $interaction_type);
        }
    }

    return @interactions;
}

sub _import_drug {
    my $self = shift;
    my $interaction = shift;
    my $citation = shift;
    my $drug_name = $self->_create_drug_name_report($interaction->{drug_id}, $citation, 'TTD_drug_id', '');

    my $primary_drug_name = $self->_create_drug_alternate_name_report($drug_name, $interaction->{drug_name}, 'TTD_primary_drug_name', '');

    my @drug_synonyms = split("; ", $interaction->{drug_synonyms});
    for my $drug_synonym (@drug_synonyms){
        next if $drug_synonym eq 'na';
        my $synonym_association = $self->_create_drug_alternate_name_report($drug_name, $drug_synonym, 'TTD_drug_synonym', '');
    }

    unless($interaction->{drug_cas_number} eq 'na'){
        my $drug_name_cas_number = $self->_create_drug_alternate_name_report($drug_name, $interaction->{drug_cas_number}, 'cas_number', '');
    }

    unless($interaction->{drug_pubchem_cid} eq 'na'){
        my $drug_name_pubchem_cid = $self->_create_drug_alternate_name_report($drug_name, $interaction->{drug_pubchem_cid}, 'pubchem_cid', '');
    }

    unless($interaction->{drug_pubchem_sid} eq 'na'){
        my $drug_name_pubchem_sid = $self->_create_drug_alternate_name_report($drug_name, $interaction->{drug_pubchem_sid}, 'pubchem_sid', '');
    }

    return $drug_name;
}

sub _import_gene {
    my $self = shift;
    my $interaction = shift;
    my $citation = shift;
    my $gene_name = $self->_create_gene_name_report($interaction->{target_id}, $citation, 'TTD_partner_id', '');

    my $gene_name_association = $self->_create_gene_alternate_name_report($gene_name, $interaction->{target_name}, 'TTD_gene_symbol', '');

    my @target_synonyms = split(";", $interaction->{target_synonyms});
    for my $target_synonym (@target_synonyms){
        next if $target_synonym eq 'na';
        my $gene_synonym = $self->_create_gene_alternate_name_report($gene_name, $target_synonym, 'TTD_alternate_gene_name', '');
    }

    my $uniprot_association = $self->_create_gene_alternate_name_report($gene_name, $interaction->{target_uniprot_id}, 'uniprot_id', '');

    return $gene_name;
}

sub input_to_tsv {
    my $self = shift;
    my $targets_raw_url = $self->targets_raw_url;
    my $drugs_crossmatching_url = $self->drugs_crossmatching_url;
    my $drug_synonyms_url = $self->drug_synonyms_url;

    #Create interactions output file
    my $interactions_outfile = $self->interactions_outfile;
    my $interactions_fh = IO::File->new($interactions_outfile, 'w');
    my $interactions_header = join("\t", 'drug_id', 'drug_name', 'drug_synonyms', 'drug_cas_number', 'drug_pubchem_cid', 'drug_pubchem_sid', 'target_id', 'target_name', 'target_synonyms', 'target_uniprot_id', 'interaction_type');
    $interactions_fh->print($interactions_header, "\n");

    #Get the data in order
    my $targets_path = $self->download_file($targets_raw_url);
    my $crossmatch_path = $self->download_file($drugs_crossmatching_url);
    my $synonyms_path = $self->download_file($drug_synonyms_url);

    my ($targets, $version) = $self->_parse_targets_file($targets_path);
    $self->version($version) if $version;
    my $drugs = $self->_parse_crossmatch_file($crossmatch_path);
    $self->_parse_synonyms_file($synonyms_path, $drugs);

    #Write data to the file
    for my $target_id (keys %{$targets}){
            #Target Uniprot Id
            my $target_uniprot_id =  pop @{$targets->{$target_id}{'UniProt ID'}};
            $target_uniprot_id = 'na' unless $target_uniprot_id;

            #Target Name
            my $target_name = pop @{$targets->{$target_id}{'Name'}};

            #Target Synonyms
            my $target_synonyms;
            $target_synonyms = join(";", @{$targets->{$target_id}{'Synonyms'}}) if $targets->{$target_id}{'Synonyms'};
            $target_synonyms = 'na' unless $target_synonyms;

        for my $drug_id (keys %{$targets->{$target_id}{'drugs'}}){
            #Drug Name
            my (undef, $drug_name) = split("\t", $targets->{$target_id}{'drugs'}{$drug_id});

            #CAS Number
            my $drug_cas_number = $drugs->{$drug_id}{'CAS Number'};
            if($drug_cas_number){
                $drug_cas_number =~ s/CAS //;
            }else{
                $drug_cas_number = 'na';
            }

            #PubChem CID
            my $drug_pubchem_cid = $drugs->{$drug_id}{'PubChem CID'};
            if($drug_pubchem_cid){
                $drug_pubchem_cid =~ s/CID //;
            }else{
                $drug_pubchem_cid = 'na';
            }

            #PubChem SID
            my $drug_pubchem_sid = $drugs->{$drug_id}{'PubChem SID'};
            if($drug_pubchem_sid){
                $drug_pubchem_sid =~ s/SID //;
            }else{
                $drug_pubchem_sid = 'na';
            }

            #Drug Synonyms
            my $drug_synonyms = $drugs->{$drug_id}{'synonyms'};
            $drug_synonyms = "na" unless $drug_synonyms;

            #Interaction Type
            my $interaction_type = $self->_determine_interaction_type($targets->{$target_id}, $drug_name, $drug_id);
            $interaction_type = 'na' unless $interaction_type;

            $interactions_fh->print(join("\t", $drug_id, $drug_name, $drug_synonyms, $drug_cas_number, $drug_pubchem_cid, $drug_pubchem_sid, $target_id, $target_name, $target_synonyms, $target_uniprot_id, $interaction_type), "\n");
        }
    }

    $interactions_fh->close;
    return 1;
}

sub preload_objects {
    my $self = shift;
    my $source_db_name = 'TTD';
    my $source_db_version = $self->version;

    #Let's preload anything for this database name and version so that we can avoid death by 1000 queries
    my @gene_names = Genome::DruggableGene::GeneNameReport->get(source_db_name => $source_db_name, source_db_version => $source_db_version);
    for my $gene_name (@gene_names){
        $gene_name->gene_alt_names;
        $gene_name->gene_categories;
    }
    my @drug_names = Genome::DruggableGene::DrugNameReport->get(source_db_name => $source_db_name, source_db_version => $source_db_version);
    for my $drug_name (@drug_names){
        $drug_name->drug_alt_names;
        $drug_name->drug_categories;
    }
    my @gene_ids = map($_->id, @gene_names);
    my @interactions = Genome::DruggableGene::DrugGeneInteractionReport->get(gene_id => \@gene_ids);
    for my $interaction (@interactions){
        $interaction->interaction_attributes;
    }

    return 1;
}

sub download_file {
    my $self = shift;
    my $url = shift;
    my ($fh, $path) = Genome::Sys->create_temp_file();
    $fh->close;
    my $wget_cmd = "wget $url -O $path";
    my $retval = Genome::Sys->shellcmd(cmd=>$wget_cmd);

    unless ($retval == 1){
      self->error_message('Failed to wget the specified URL');
      return;
    }

    return $path;
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
