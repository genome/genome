package Genome::Db::Ensembl::AnnotationStructures;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use IO::File;
use Carp;
use Sys::Hostname;

class Genome::Db::Ensembl::AnnotationStructures {
    is  => 'Genome::SoftwareResult::Stageable',
    has => [
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'reference_build_id',
        },
        log_file => {
            is => 'Text',
            is_optional  => 1,
            is_transient => 1,
        },
        dump_file => {
            is => 'Text',
            is_optional  => 1,
            is_transient => 1,
        },
    ],
    has_param => [
        version => {
            is  => 'Text',
            doc => "Version to use",
        },
        species => {
            is => 'Text',
            doc => 'Species of annotation to import (mouse, human currently supported)',
            valid_values => [qw(mouse human)],
        },
        data_set => {
            is => 'Text',
            doc => 'Ensembl data set to import',
            default => 'Core',
        },
        software_version => {
            is => 'Text',
            doc => 'This should be incremented when a software change is made that would change the annotation structures',
        },
    ],
    has_input =>  [
        reference_build_id => {
            is => 'Text',
        }
    ],
    has_optional_transient => [
        _user_test_name => {
            is => 'Text',
            doc => "If create is called with a test name, don't unset it at the end",
        }
    ],
};

sub create
{
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    if ($self->test_name) {
        $self->_user_test_name($self->test_name);
    }
    else {
        #This is necessary because we call commit after each chromosome,
        #so even if it crashes, the software result will stick around
        $self->set_test_name("Creation of this software result is not complete");
    }

    $self->_prepare_staging_directory;
    $self->status_message('Create AnnotationStructures');
    my $data_set = $self->data_set;
    $self->prepare_for_execution;

    my $registry = $self->connect_registry;
    my $db_connection = $self->get_db_connection($registry);
    my $ucfirst_species = ucfirst $self->species;
    my $gene_adaptor = $registry->get_adaptor( $ucfirst_species, $data_set, 'Gene' );
    my $transcript_adaptor = $registry->get_adaptor( $ucfirst_species, $data_set, 'Transcript' );
    my $slice_adaptor = $registry->get_adaptor( $ucfirst_species, $data_set, 'Slice');
    my $assembly_exception_feature_adaptor = $registry->get_adaptor( $ucfirst_species, $data_set, 'AssemblyExceptionFeature');

    # Map Ensembl and GENCODE biotypes to more general categories that our annotator can use
    my %biotype_class;
    map{$biotype_class{lc($_)} = 'ncrna'} qw( 3prime_overlapping_ncrna Mt_rRNA Mt_tRNA antisense lincRNA miRNA misc_RNA ncRNA ncrna_host non_coding processed_transcript rRNA sense_intronic sense_overlapping snRNA snoRNA tRNA );
    map{$biotype_class{lc($_)} = 'ncrna_pseudogene'} qw( Mt_tRNA_pseudogene miRNA_pseudogene misc_RNA_pseudogene rRNA_pseudogene scRNA_pseudogene snRNA_pseudogene snoRNA_pseudogene tRNA_pseudogene );
    map{$biotype_class{lc($_)} = 'coding'} qw( LRG_gene IG_C_gene IG_D_gene IG_J_gene IG_V_gene TR_C_gene TR_D_gene TR_J_gene TR_V_gene protein_coding ambiguous_orf disrupted_domain non_stop_decay nonsense_mediated_decay retained_intron TEC );
    map{$biotype_class{lc($_)} = 'coding_pseudogene'} qw( retrotransposed TR_J_pseudogene TR_V_pseudogene IG_C_pseudogene IG_J_pseudogene IG_V_pseudogene polymorphic_pseudogene processed_pseudogene pseudogene transcribed_processed_pseudogene transcribed_unprocessed_pseudogene unitary_pseudogene unprocessed_pseudogene );

    my @slices = @{ $slice_adaptor->fetch_all('toplevel', undef, 1, 1, 1) };

    my $idx = 0;
    my $egi_id = 1;    # starting point for external_gene_id...
    my $tss_id = 1;    # starting point for transcript sub struct ids...

    my $count = 0;
    my $pfcount = 0;
    #species, source and version are id properties on all items
    my $source = 'ensembl';
    my $version = $self->version;
    my $species = $self->species;

    #for logging purposes
    my @transcripts;
    my @sub_structures;
    my @proteins;
    my @genes;

    foreach my $slice ( @slices ) {
        my $chromosome = $slice->seq_region_name();
        my $bases_file = $self->reference_build->get_bases_file($chromosome);
        unless (-e $bases_file) {
            $self->status_message("Transcript structures from chromosome $chromosome will not be imported because it doesn\'t exist in the reference sequence.");
            next;
        }
        $self->status_message("Importing transcripts for $chromosome");
        my @ensembl_transcripts = @{ $transcript_adaptor->fetch_all_by_Slice($slice)};
        $self->status_message("Importing ".scalar @ensembl_transcripts." transcripts\n");

        #We are included duplicate regions (like PARs and HAPs) in order
        #to get the coordinate system right.  However, we don't actually
        #want to import duplicate transcripts.
        #Keep track of which regions are duplicates and don't import
        #those transcripts.
        my @assembly_exceptions = @{ $assembly_exception_feature_adaptor->fetch_all_by_Slice($slice)};
        #Ensembl api returns assembly exception two assembly exception features
        #for each assembly exception, one going each way.  We only want to
        #exclude transcripts from the non-REF copy
        @assembly_exceptions = grep {!($_->type =~ /REF/)} @assembly_exceptions;

        TRANSCRIPT: foreach my $ensembl_transcript (@ensembl_transcripts) {
            foreach my $assembly_exception (@assembly_exceptions) {
                if ($assembly_exception->overlaps($ensembl_transcript)) {
                    my $print_name = $ensembl_transcript->stable_id;
                    $self->warning_message("Skipping transcript $print_name because it is on an assembly_exception");
                    next TRANSCRIPT;
                }
            }

            $count++;

            #The Ensembl/GENCODE biotype is used to determine RNA sub_structures and pseudogene status
            my $biotype = $ensembl_transcript->biotype();
            unless (defined $biotype_class{lc($biotype)}) {
                $self->error_message("Novel biotype \"" . $biotype . "\" seen in DB! Please update the \%biotype_class hash.");
                die;
            }
            $biotype = lc( $biotype );

            my $ensembl_gene  = $gene_adaptor->fetch_by_transcript_id( $ensembl_transcript->dbID );
            my $ensembl_gene_id = $ensembl_gene->dbID;

            my $transcript_start = $ensembl_transcript->start;
            my $transcript_stop  = $ensembl_transcript->end;

            next unless defined $transcript_start and defined $transcript_stop;

            my $strand = $ensembl_transcript->strand;
            if ( $strand == 1 ) {
                $strand = "+1";
            }
            elsif ($strand == -1) {
                $strand = "-1";
            }
            else {
                $self->warning_message("Invalid strand $strand, skipping!");
            }

            my $hugo_gene_name = undef;
            my $external_db;
            if ($species eq 'human') {
                $external_db = 'HGNC';
            }
            elsif($species eq 'mouse') {
                $external_db = 'MGI';
            }
            if ($ensembl_gene->external_db =~ /$external_db/){
                $hugo_gene_name = $ensembl_gene->external_name;
            }

            my $entrez_id = undef;
            my $entrez_genes = $ensembl_gene->get_all_DBEntries('EntrezGene');
            if ( defined(@$entrez_genes)) {
                $entrez_id = @$entrez_genes[0]->primary_id;
            }

            #gene cols: gene_id hugo_gene_name strand
            my $gene;
            my $gene_meta = Genome::Gene->__meta__;
            my $composite_gene_id = $gene_meta->resolve_composite_id_from_ordered_values($ensembl_gene_id, $species,$source,$version);
            $gene = Genome::Gene->get(id => $composite_gene_id, data_directory => $self->temp_staging_directory, reference_build_id => $self->reference_build->id);
            unless ($gene) {
                $gene = Genome::Gene->create(
                    gene_id => $ensembl_gene_id,
                    hugo_gene_name => $hugo_gene_name,
                    strand => $strand,
                    data_directory => $self->temp_staging_directory,
                    species => $species,
                    source => $source,
                    version => $version,
                    reference_build_id => $self->reference_build->id,
                );
                push @genes, $gene;#logging

                my %external_gene_ids = $self->get_external_gene_ids($ensembl_gene);
                if ( defined($hugo_gene_name) ) {
                    $external_gene_ids{hugo_symbol} = $hugo_gene_name;
                }

                if ( defined($entrez_id) ) {
                    $external_gene_ids{entrez} = $entrez_id;
                }
                $external_gene_ids{ensembl} = $ensembl_gene->stable_id;

                #external_gene_id columns
                foreach my $type ( sort keys %external_gene_ids ) {
                    my $external_gene_id = Genome::ExternalGeneId->create(
                        egi_id => $egi_id,
                        gene_id => $gene->id,
                        id_type => $type,
                        id_value => $external_gene_ids{$type},
                        data_directory => $self->temp_staging_directory,
                        species => $species,
                        source => $source,
                        version => $version,
                        reference_build_id => $self->reference_build->id,
                    );

                    $egi_id++;
                }
                my $external_gene_id = Genome::ExternalGeneId->create(
                    egi_id => $egi_id,
                    gene_id => $gene->id,
                    id_type => "ensembl_default_external_name",
                    id_value => $ensembl_gene->external_name,
                    data_directory => $self->temp_staging_directory,
                    species => $species,
                    source => $source,
                    version => $version,
                    reference_build_id => $self->reference_build->id,
                );

                $egi_id++;
                $external_gene_id = Genome::ExternalGeneId->create(
                    egi_id => $egi_id,
                    gene_id => $gene->id,
                    id_type => "ensembl_default_external_name_db",
                    id_value => $ensembl_gene->external_db,
                    data_directory => $self->temp_staging_directory,
                    species => $species,
                    source => $source,
                    version => $version,
                    reference_build_id => $self->reference_build->id,
                );

                $egi_id++;
            }

            #Transcript cols: transcript_id gene_id transcript_start transcript_stop transcript_name source transcript_status strand chrom_name

            my $transcript = Genome::Transcript->create(
                reference_build_id => $self->reference_build_id,
                transcript_id => $ensembl_transcript->dbID,
                gene_id => $gene->id,
                gene_name => $gene->name,
                transcript_start => $transcript_start,
                transcript_stop => $transcript_stop,
                transcript_name => $ensembl_transcript->stable_id,
                transcript_status => lc( $ensembl_transcript->status ),   #TODO valid statuses (unknown, known, novel) #TODO verify substructures and change status if necessary
                strand => $strand,
                chrom_name => $chromosome,
                data_directory => $self->temp_staging_directory,
                species => $species,
                source => $source,
                version => $version,
            );
            push @transcripts, $transcript; #logging

            #sub structures
            my @ensembl_exons = @{ $ensembl_transcript->get_all_Exons() };
            my @utr_exons;
            my @cds_exons;
            my @rna;

            foreach my $exon ( @ensembl_exons ) {
                #Ensembl exons are combined CDS and UTR regions, we need to create both utr and cds exons.
                #From these, flank and intron will be created after instantianting these substructures.
                #Phase and ordinal will be set after instantiating these substructures
                my $coding_region_start = $exon->coding_region_start($ensembl_transcript);
                my $coding_region_stop = $exon->coding_region_end($ensembl_transcript);
                my $start = $exon->start;
                my $stop = $exon->end;
                my $exon_sequence = $exon->seq->seq;

                if (defined $coding_region_start) {
                    #There is a coding section in this exon
                    unless (defined $coding_region_stop) {
                        $self->error_message("ensembl exon has a coding_region_start defined, but not a coding_region_end!". Dumper $exon);
                        die;
                    }

                    if ($coding_region_start > $start) {
                        #there is a utr exon between the start of the transcript and the coding region
                        #create utr_exon
                        my $utr_sequence;
                        if ( $transcript->strand eq '+1' ) {
                            $utr_sequence = substr( $exon_sequence, 0, $coding_region_start - $start );
                        }
                        else {
                            #sequence is returned stranded, so we need the seq from the end
                            $utr_sequence = substr( $exon_sequence, 0 - ( $coding_region_start - $start ))
                        }

                        my $utr_stop = $coding_region_start - 1;
                        my $utr_exon = $self->create_transcript_structure(
                            transcript => $transcript,
                            chrom_name => $transcript->chrom_name,
                            transcript_structure_id => $tss_id,
                            transcript_id => $transcript->id,
                            structure_type => 'utr_exon',
                            structure_start => $start,
                            structure_stop => $utr_stop,
                            nucleotide_seq => $utr_sequence,
                            data_directory => $self->temp_staging_directory,
                            species => $species,
                            source => $source,
                            version => $version,
                        );
                        $tss_id++;
                        push @utr_exons, $utr_exon;
                        push @sub_structures, $utr_exon; #logging
                    }

                    #create cds_exon (we do a little extra arithmetic here if the whole exon is coding, cleaner this way but could add an alternative block if coding_region_start == start and coding_region_stop == stop)
                    my $cds_sequence;
                    if ( $transcript->strand eq '+1' ) {
                        #grab sequence from start of coding region for length of coding region
                        $cds_sequence = substr(
                            $exon_sequence,
                            $coding_region_start - $start,
                            $coding_region_stop - $coding_region_start + 1 );
                    }
                    else {
                        #otherwise grab starting at the index of the distance from the coding_region stop to the stop for the length of the coding region
                        $cds_sequence
                        = substr( $exon_sequence, $stop - $coding_region_stop, $coding_region_stop - $coding_region_start + 1 );
                    }

                    my $cds_exon = $self->create_transcript_structure(
                        transcript => $transcript,
                        chrom_name => $transcript->chrom_name,
                        transcript_structure_id => $tss_id,
                        transcript_id => $transcript->id,
                        structure_type => 'cds_exon',
                        structure_start => $coding_region_start,
                        structure_stop => $coding_region_stop,
                        nucleotide_seq => $cds_sequence,
                        data_directory => $self->temp_staging_directory,
                        species => $species,
                        source => $source,
                        version => $version,
                    );

                    $tss_id++;
                    push @cds_exons, $cds_exon;
                    push @sub_structures, $cds_exon; #logging

                    if ($stop > $coding_region_stop) {
                        #there is a utr exon after the coding region
                        #create utr_exon
                        my $utr_sequence;
                        if ( $transcript->strand eq '+1' ) {
                            $utr_sequence = substr( $exon_sequence, 0 - ( $stop - $coding_region_stop ));
                        }
                        else {
                            $utr_sequence = substr( $exon_sequence, 0, $stop - $coding_region_stop )
                        }

                        my $utr_start = $coding_region_stop + 1;

                        my $utr_exon = $self->create_transcript_structure(
                            transcript => $transcript,
                            chrom_name => $transcript->chrom_name,
                            transcript_structure_id => $tss_id,
                            transcript_id => $transcript->id,
                            structure_type => 'utr_exon',
                            structure_start => $utr_start,
                            structure_stop => $stop,
                            nucleotide_seq => $utr_sequence,
                            data_directory => $self->temp_staging_directory,
                            species => $species,
                            source => $source,
                            version => $version,
                        );
                        $tss_id++;
                        push @utr_exons, $utr_exon;
                        push @sub_structures, $utr_exon; #logging

                    }
                }
                elsif(defined $coding_region_stop) {
                    $self->error_message("ensembl exon has a coding_region_end, but not a coding_region_start!". Dumper $exon);
                    die;
                }
                else {
                    #If we reached here, then this entire exon is either a UTR or non-coding RNA
                    my $structure_type = 'utr_exon';

                    #Check if this exon belongs to a non-coding transcript. Include ncRNA pseudogenes
                    if ($biotype_class{$biotype} eq 'ncrna' or $biotype_class{$biotype} eq 'ncrna_pseudogene') {
                        $structure_type = 'rna';
                    }

                    my $structure = $self->create_transcript_structure(
                        transcript => $transcript,
                        chrom_name => $transcript->chrom_name,
                        transcript_structure_id => $tss_id,
                        transcript_id => $transcript->id,
                        structure_type => $structure_type,
                        structure_start => $start,
                        structure_stop => $stop,
                        nucleotide_seq => $exon_sequence,
                        data_directory => $self->temp_staging_directory,
                        species => $species,
                        source => $source,
                        version => $version,
                    );

                    push @rna, $structure if $structure_type eq 'rna';
                    push @utr_exons, $structure unless $structure_type eq 'rna';
                    $tss_id++;
                    push @sub_structures, $structure;
                }
            }
            if (@utr_exons > 0 or @cds_exons > 0) {
                $self->assign_ordinality_to_exons( $transcript->strand, [@utr_exons, @cds_exons] );
            }
            if (@cds_exons > 0) {
                $self->assign_phase( \@cds_exons );
            }

            #create flanks and intron
            my @flanks_and_introns = $self->create_flanking_sub_structures_and_introns($transcript, \$tss_id, [@cds_exons, @utr_exons, @rna]);

            my $protein;
            my $translation = $ensembl_transcript->translation();
            if ( defined($translation) ) {
                $protein = Genome::Protein->create(
                    protein_id => $translation->dbID,
                    transcript_id => $transcript->id,
                    protein_name => $translation->stable_id,
                    amino_acid_seq => $ensembl_transcript->translate->seq,
                    data_directory => $self->temp_staging_directory,
                    species => $species,
                    source => $source,
                    version => $version,
                    reference_build_id => $self->reference_build->id,
                );
                push @proteins, $protein;

                my $protein_features = $translation->get_all_ProteinFeatures;

                foreach my $pf (@$protein_features) {
                    $pfcount++;
                    my $logic_name = $pf->analysis->logic_name();
                    my $pf_start = $pf->start;
                    my $pf_end = $pf->end;
                    my $interpro_ac = $pf->interpro_ac;
                    my $idesc = $pf->idesc;
                    my $program = $pf->analysis->db;

                    if (!$interpro_ac) {
                        if ($pf->hseqname =~ /SSF/) {
                            #try prepending the superfamily prefix and querying directly
                            my $new_id = $pf->hseqname;
                            $new_id =~ s/SSF//;
                            my $query = "select i.interpro_ac,x.display_label from interpro i join xref x on i.interpro_ac = x.dbprimary_acc where i. id=$new_id;";
                            my $sth = $db_connection->prepare($query);
                            my $exec_result = $sth->execute;
                            if (!$exec_result) {
                                #query failed
                                $self->error_message("Failed to execute sql query for interpro: $query");
                                return;
                            }
                            if ($exec_result == 0) {
                                $self->warning_message("Transcript had protein domain, but we couldn't get the domain name out of ensembl.  Protein feature hit name: ".$pf->hseqname);
                                next;
                            }
                            my $result = $sth->fetchrow_hashref;
                            $interpro_ac = $result->{interpro_ac};

                            $idesc = $result->{display_label};
                        }
                        else {
                            $self->warning_message("Transcript had protein domain, but ensembl didn't have a domain name.  Protein feature hit name: ".$pf->hseqname);
                                next;
                        }
                    }
                    unless ($program) {$program = 'na'};

                    my $interpro_result = Genome::InterproResult->create(
                        interpro_id => $pfcount,
                        chrom_name => $chromosome,
                        transcript_name => $transcript->transcript_name,
                        data_directory => $self->temp_staging_directory,
                        start => $pf_start,
                        stop => $pf_end,
                        rid => 0, #Copied directly from mg-load-ipro
                        setid => 'na',
                        parent_id => 'na',
                        name => $logic_name."_".$idesc,
                        interpro_note => "Interpro ".$interpro_ac." ".$idesc,
                    );
                }
            }

            if ($transcript->cds_full_nucleotide_sequence) {
                my $transcript_seq = Genome::TranscriptCodingSequence->create(
                    transcript_id => $transcript->id,
                    sequence => $transcript->cds_full_nucleotide_sequence,
                    data_directory => $transcript->data_directory,
                );
            }

            # Assign various fields to the transcript
            my %transcript_info;
            $transcript_info{pseudogene} = 0;
            if ( $biotype_class{$biotype} eq 'ncrna_pseudogene' or $biotype_class{$biotype} eq 'coding_pseudogene' ) {
                $transcript_info{pseudogene} = 1;
            }
            $self->calculate_transcript_info($transcript, \%transcript_info);

            my @structures = $transcript->sub_structures;
            foreach my $structure (@structures) {
                $self->_update_transcript_info($structure, $transcript);
            }
        }

        $self->write_log_entry($count, \@transcripts, \@sub_structures, \@genes, \@proteins);

        #$self->dump_sub_structures(0); #arg added for pre/post commit notation

        $self->status_message( "committing...($count)" );
        UR::Context->commit;
        $self->status_message( "finished commit!\n" );

        #$self->dump_sub_structures(1);

        Genome::Gene->unload;
        Genome::Transcript->unload;
        Genome::ExternalGeneId->unload;
        Genome::TranscriptStructure->unload;
        Genome::Protein->unload;
        Genome::TranscriptCodingSequence->unload;
        Genome::InterproResult->unload;

        #reset logging arrays
        @transcripts = ();
        @genes = ();
        @proteins = ();
        @sub_structures = ();

        #exit; #uncomment for testing
    }
    $self->status_message( "Create AnnotationStructures done" );

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    unless ($self->_user_test_name) {
        $self->remove_test_name();
    }
    UR::Context->commit;

    return $self;
}

sub result_paths {
    return ("external_gene_ids.csv", "genes", "interpro_results", "proteins", "substructures", "transcript_coding_sequences.csv", "transcripts.csv");
}

#TODO: make this connect to a public ensembl DB if the environment variables      aren't set
sub get_ensembl_info {
    my $self = shift;

    my $host = defined $ENV{GENOME_DB_ENSEMBL_HOST} ?                             $ENV{GENOME_DB_ENSEMBL_HOST} : 'mysql1';
    my $user = defined $ENV{GENOME_DB_ENSEMBL_USER} ?                             $ENV{GENOME_DB_ENSEMBL_USER} : 'mse';
    my $password = defined $ENV{GENOME_DB_ENSEMBL_PASSWORD} ?                     $ENV{GENOME_DB_ENSEMBL_PASSWORD} : undef;

    return ($host, $user, $password);
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("ensembl_annotation_structures-%s-%s-%s-%s",           $hostname, $user, $$, $self->id);
    my $directory = join('/', 'build_merged_alignments',$self->id,$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    return 'info_genome_models';
}

sub _update_transcript_info {
    my $self = shift;
    my $structure = shift;
    my $transcript = shift;

    $structure->transcript_gene_name($transcript->gene_name);
    $structure->transcript_transcript_error($transcript->transcript_error);
    $structure->transcript_coding_region_start($transcript->coding_region_start);
    $structure->transcript_coding_region_stop($transcript->coding_region_stop);
    $structure->transcript_amino_acid_length($transcript->amino_acid_length);

    return 1;
}

sub connect_registry {
    my $self = shift;

    my $lib = "use Bio::EnsEMBL::Registry;\nuse Bio::EnsEMBL::DBSQL::DBAdaptor;";
    eval $lib;

    if ($@) {
        $self->error_message("not able to load the ensembl modules");
        croak();
    }

    my ($host, $user, $pass) = $self->get_ensembl_info($self->version);

    my $registry = 'Bio::EnsEMBL::Registry';
    my $ens_host = $host;
    my $ens_user = $user;

    $registry->load_registry_from_db(
        -host => $ens_host,
        -user => $ens_user,
    );
    return $registry;
}

sub get_db_connection {
    my $self = shift;
    my $registry = shift;
    my $ucfirst_species = ucfirst $self->species;
    my $db_adaptor = $registry->get_DBAdaptor($ucfirst_species, "core");
    my $db_connection = $db_adaptor->dbc;
    return $db_connection;
}

sub ordcount {
    my $self = shift;
    my $ord  = shift;
    my $type = shift;
    if ( !defined( $ord->{$type} )) {
        $ord->{$type} = 1;
    }
    else {
        $ord->{$type}++;
    }
    return $ord->{$type};
}

sub get_external_gene_ids {
    my $self = shift;
    my $gene = shift;
    my %external_ids;
    my @entries = @{$gene->get_all_DBEntries()};
    my @dbswanted = qw/ UCSC EntrezGene OTTT CCDS Vega_gene /;
    if ($self->species eq 'human') {
        unshift @dbswanted, (qw/HGNC HGNC_automatic_gene/);
    }
    elsif ($self->species eq 'mouse') {
        unshift @dbswanted, (qw/MGI MGI_automatic_gene/);
    }
    my %dbs = map { $_ => 1 } @dbswanted;
    foreach my $entry (@entries) {
        my $dbname = $entry->dbname();
        next unless exists($dbs{$dbname});
        my $dbvalue = $entry->display_id();
        $external_ids{$dbname} = $dbvalue;
    }
    return %external_ids;
}

sub prepare_for_execution {
    my $self = shift;
    if (defined $self->log_file and -e $self->log_file) {
        unlink $self->log_file;
    }
    if (defined $self->dump_file and -e $self->dump_file) {
        unlink $self->dump_file;
    }

    Genome::Sys->create_directory($self->temp_staging_directory."/genes");
    Genome::Sys->create_directory($self->temp_staging_directory."/substructures");
    Genome::Sys->create_directory($self->temp_staging_directory."/proteins");
    Genome::Sys->create_directory($self->temp_staging_directory."/interpro_results");

    $self->{cumulative_transcripts} = 0;
    $self->{cumulative_sub_structures} = 0;
    $self->{cumulative_genes} = 0;
    $self->{cumulative_proteins} = 0;

    return 1;
}

sub assign_ordinality {
    my ($self, $strand, $ss_array) = @_;
    my @a;
    if ($strand eq '+1'){
        @a = sort {$a->structure_start <=> $b->structure_start} @$ss_array;
    }
    else {
        @a =  sort {$b->structure_start <=> $a->structure_start} @$ss_array;
    }
    my $ord = 1;
    for my $ss (@a){
        $ss->ordinal($ord);
        $ord++;
    }
    return 1;
}

sub assign_ordinality_to_exons {
    my ($self, $strand, $ss_array) = @_;
    my @a;
    if ($strand eq '+1') {
        @a = sort {$a->structure_start <=> $b->structure_start} @$ss_array;
    }
    else {
        @a =  sort {$b->structure_start <=> $a->structure_start} @$ss_array;
    }
    my $ord = 1;
    my $last_exon;
    for my $ss (@a) {
        unless ($last_exon) {
            $ss->ordinal($ord);
            $last_exon = $ss;
            next;
        }
        $ord++ unless $self->sub_structures_are_contiguous($strand, $last_exon, $ss);
        $ss->ordinal($ord);
        $last_exon = $ss;
    }
    return 1;
}

sub sub_structures_are_contiguous {
    my ($self, $strand, $a, $b) = @_;
    if ($strand eq '+1') {
        if ($a->structure_stop >= $b->structure_start) {
            $self->warning_message('exons overlap!');
            return undef;
        }
        if ($a->structure_stop+1 == $b->structure_start) {
            return 1;
        }

    }
    elsif($strand eq '-1') {
        if ($b->structure_stop >= $a->structure_start) {
            $self->warning_message('exons overlap!');
            return undef;
        }
        if ($b->structure_stop+1 == $a->structure_start) {
            return 1;
        }
    }
}

sub assign_phase {
    my ($self, $cds_array) = @_;
    my @a = sort {$a->ordinal <=> $b->ordinal} @$cds_array;
    my $phase = 0;

    my $previous_exon = shift @a;
    $previous_exon->phase($phase);

    for my $cds_exon (@a) {
        $phase = ( $phase + $previous_exon->length ) % 3;
        $cds_exon->phase($phase);
        $previous_exon = $cds_exon;
    }

    return 1;
}

# Calls several other methods in this module that calculate particular fields for the transcript
sub calculate_transcript_info {
    my ($self, $transcript, $transcript_info) = @_;

    my @methods = qw/
    phase_nucleotides
    coding_bases_before_and_after_substructures
    coding_exons_before_and_after_substructures
    calculate_transcript_coding_region
    calculate_transcript_amino_acid_length
    determine_error_status /;

    for my $method (@methods) {
        my $rv = $self->$method($transcript);
        confess "Problem executing $method during annotation import!" unless $rv;
    }

    return 1 unless defined $transcript_info;

    # Handle various params passed in transcript info hash
    # Allows data available to the importer to be used to calculate and store transcript info
    if (exists $transcript_info->{pseudogene}) {
        if ($transcript_info->{pseudogene} and $transcript->internal_stop_codon == 1) {
            my $error = $transcript->transcript_error;
            $error = "" if $error eq 'no_errors';
            $transcript->transcript_error($error . ":pseudogene") unless $error eq "";
            $transcript->transcript_error("pseudogene") if $error eq "";
        }
    }

    return 1;
}

# Exons that use nucleotides from previous/next exons to complete a codon     (ie, have nonzero phase)
# have those nucleotides stored so they don't have to pull the full sequence  from neighbors
sub phase_nucleotides {
    my ($self, $transcript) = @_;
    my @ss = $transcript->ordered_sub_structures;

    my $previous_exon;
    while (my $s = shift @ss) {
        unless ($s->{structure_type} eq 'cds_exon') {
            $s->{phase_bases_before} = 'NULL';
            $s->{phase_bases_after} = 'NULL';
            next;
        }

        if ($s->phase != 0) {
            unless (defined $previous_exon) {
                $self->error_message("Phase of exon " . $s->{transcript_structure_id} . " on transcript " .
                    $transcript->{transcript_name} . " indicates that bases from previous exon needed to " .
                    "complete first codon, but there is no previous exon!");
                confess;
            }

            $s->{phase_bases_before} = substr($previous_exon->nucleotide_seq, -1 * $s->phase);
        }
        else {
            $s->{phase_bases_before} = 'NULL';
        }

        my $next_exon;
        for my $next (@ss) { $next_exon = $next and last if $next->{structure_type} eq 'cds_exon' }

        if ($next_exon and $next_exon->phase != 0) {
            $s->{phase_bases_after} = substr($next_exon->nucleotide_seq, 0, 3 - $next_exon->phase);
        }
        else {
            $s->{phase_bases_after} = 'NULL';
        }
        $previous_exon = $s;
        # Construct sequence and make sure it contains no partial codons (should be divisible by 3)
        my $seq;
        unless ($s->{phase_bases_before} eq 'NULL') {
            $seq .= $s->{phase_bases_before};
        }
        $seq .= $s->nucleotide_seq;
        unless ($s->{phase_bases_after} eq 'NULL') {
            $seq .= $s->{phase_bases_after};
        }

        unless ((length $seq) % 3 == 0) {
            $self->warning_message("Sequence constructed using phase bases for exon " . $s->{transcript_structure_id} .
                " from transcript " . $transcript->{transcript_name} . " contains an incomplete codon!");
        }
    }

    return 1;
}

# All substructures need to know how many coding basepairs exist both before  and after them in the transcript
# to ease calculation of the coding sequence position during annotation
sub coding_bases_before_and_after_substructures {
    my ($self, $transcript) = @_;

    my @ss = $transcript->ordered_sub_structures;
    for my $s (@ss) {
        my $coding_bases_before = $transcript->length_of_cds_exons_before_structure_at_position(
            $s->structure_start,
        );
        unless (defined $coding_bases_before) {
            $self->error_message("Could not calculate previous coding basepairs for substructure "
                . $s->transcript_structure_id . " on transcript " . $transcript->{transcript_name});
            confess;
        }

        my $coding_bases_after = $transcript->length_of_cds_exons_after_structure_at_position(
            $s->structure_start,
        );
        unless (defined $coding_bases_after) {
            $self->error_message("Could not calculate following coding basepairs for substructure " .
                $s->transcript_structure_id . " on transcript " . $transcript->{transcript_name});
            confess;
        }

        $s->{coding_bases_before} = $coding_bases_before;
        $s->{coding_bases_after} = $coding_bases_after;
    }

    return 1;
}

# All substructures need to know how many coding exons exists before and      after them in the transcript
# to further ease annotation
sub coding_exons_before_and_after_substructures {
    my ($self, $transcript) = @_;
    my @ss = $transcript->ordered_sub_structures;
    my $exons_passed = 0;
    my $unreached_exons = scalar $transcript->cds_exons;

    for my $s (@ss) {
        $unreached_exons-- if $s->{structure_type} eq 'cds_exon';
        $s->{cds_exons_before} = $exons_passed;
        $s->{cds_exons_after} = $unreached_exons;
        $exons_passed++ if $s->{structure_type} eq 'cds_exon';
    }

    return 1;
}

# Assigns the start and stop position of the transcript's coding region
sub calculate_transcript_coding_region {
    my ($self, $transcript) = @_;
    my ($coding_region_start, $coding_region_stop) = $transcript->cds_exon_range;
    unless (defined $coding_region_start and defined $coding_region_stop) {
        $transcript->{coding_region_start} = 'NULL';
        $transcript->{coding_region_stop} = 'NULL';
        return 1;
    }

    $transcript->{coding_region_start} = $coding_region_start;
    $transcript->{coding_region_stop} = $coding_region_stop;
    return 1;
}

# Assigns the length of the amino acid to the transcript
sub calculate_transcript_amino_acid_length {
    my ($self, $transcript) = @_;
    my $amino_acid_length = 0;
    my $protein = $transcript->protein;
    $amino_acid_length = length $protein->amino_acid_seq if $protein;

    $transcript->{amino_acid_length} = $amino_acid_length;
    return 1;
}

# Set transcript error field
sub determine_error_status {
    my ($self, $transcript) = @_;
    $transcript->is_valid;
    return 1;
}

sub create_flanking_sub_structures_and_introns {
    my ($self, $transcript, $tss_id_ref, $ss_array) = @_;
    return unless @$ss_array;

    my @a = sort {$a->structure_start <=> $b->structure_start} @$ss_array;

    my $left_flank_structure_stop = $a[0]->structure_start - 1;
    my $left_flank_structure_start = $a[0]->structure_start - 50000;
    my $left_flank = $self->create_transcript_structure(
        transcript => $transcript,
        chrom_name => $transcript->chrom_name,
        transcript_structure_id => $$tss_id_ref,
        transcript_id => $transcript->id,
        structure_type => 'flank',
        structure_start => $left_flank_structure_start,
        structure_stop => $left_flank_structure_stop,
        data_directory => $self->temp_staging_directory,
        species => $self->species,
        source => $transcript->source,
        version => $self->version,

    );
    $$tss_id_ref++;

    my $right_flank_structure_start = $a[-1]->structure_stop + 1;
    my $right_flank_structure_stop = $a[-1]->structure_stop + 50000;
    my $right_flank = $self->create_transcript_structure(
        transcript => $transcript,
        chrom_name => $transcript->chrom_name,
        transcript_structure_id => $$tss_id_ref,
        transcript_id => $transcript->id,
        structure_type => 'flank',
        structure_start => $right_flank_structure_start,
        structure_stop => $right_flank_structure_stop,
        data_directory => $self->temp_staging_directory,
        species => $self->species,
        source => $transcript->source,
        version => $self->version,
    );
    $$tss_id_ref++;

    $self->assign_ordinality($transcript->strand, [$left_flank,
        $right_flank]);

    #now create introns for any gaps between exons
    my @introns;
    my $left_ss;
    for my $ss (@a) {
        unless ($left_ss) {
            $left_ss = $ss;
            next;
        }
        my $right_structure_start = $ss->structure_start;
        my $left_structure_stop = $left_ss->structure_stop;

        if ( $right_structure_start > $left_structure_stop + 1 ) {
            my $intron_start = $left_structure_stop+1;
            my $intron_stop = $right_structure_start-1;
            my $intron = $self->create_transcript_structure (
                transcript => $transcript,
                chrom_name => $transcript->chrom_name,
                transcript_structure_id => $$tss_id_ref,
                transcript_id => $transcript->id,
                structure_type => 'intron',
                structure_start => $intron_start,
                structure_stop => $intron_stop,
                data_directory => $self->temp_staging_directory,
                species => $self->species,
                source => $transcript->source,
                version => $self->version,
            );
            $$tss_id_ref++;
            push @introns, $intron
        }
        $left_ss = $ss;
    }
    $self->assign_ordinality($transcript->strand, \@introns);
    return ($left_flank, $right_flank, @introns);
}

#Convenience method to tag on all of the transcript fields
sub create_transcript_structure {
    my ($class, %params) = @_;
    my $transcript = delete $params{transcript};

    map {if ($transcript->$_) {my $param_name = 'transcript_'.$_;
    $params{$param_name} = $transcript->$_}} qw/transcript_id gene_id transcript_start transcript_stop transcript_name transcript_status strand chrom_name species source version gene_name transcript_error coding_region_start coding_region_stop amino_acid_length/;

    return Genome::TranscriptStructure->create(%params);
}

sub write_log_entry {
    my $self = shift;
    return unless defined $self->log_file;
    my ($count, $transcripts, $sub_structures, $genes, $proteins) = @_;
    my $log_fh = IO::File->new(">> ". $self->log_file);
    return unless $log_fh;
    my $total_transcripts = scalar @$transcripts;
    $self->{cumulative_transcripts} += $total_transcripts;
    my $total_sub_structures = scalar @$sub_structures;
    $self->{cumulative_sub_structures} += $total_sub_structures;
    my $total_genes = scalar @$genes;
    $self->{cumulative_genes} += $total_genes;
    my $total_proteins = scalar @$proteins;
    $self->{cumulative_proteins} += $total_proteins;
    my @transcripts_missing_sub_structures;
    my @transcripts_missing_gene;
    for my $transcript (@$transcripts) {
        my @ss = $transcript->sub_structures;
        my $gene = $transcript->gene;
        push @transcripts_missing_sub_structures, $transcript unless @ss;
        push @transcripts_missing_gene, $transcript unless $gene;
    }

    $log_fh->print("Count $count\n");
    $log_fh->print("transcripts this round $total_transcripts\n");
    $log_fh->print("sub_structures this round $total_sub_structures\n");
    $log_fh->print("genes this round $total_genes\n");
    $log_fh->print("proteins this round $total_proteins\n");
    $log_fh->print("Cumulative: ".$self->{cumulative_transcripts}."           transcripts, ".$self->{cumulative_sub_structures}." ss, ".$self->{cumulative_genes}." genes, ".$self->{cumulative_proteins}." proteins\n");
    $log_fh->print("There were ".scalar @transcripts_missing_sub_structures." transcripts missing sub_structures: ".join(" ", map {$_->transcript_name} @transcripts_missing_sub_structures)."\n");
    $log_fh->print("There were ".scalar @transcripts_missing_gene."           transcripts missing a gene ".join(" ", map {$_->transcript_name} @transcripts_missing_gene)."\n");
    $log_fh->print("##########################################\n");
}

sub dump_sub_structures {
    my $self = shift;
    my ($committed) = @_; #indicates if we are dumping status pre or post     commit
    my $dump_fh = IO::File->new(">> ". $self->dump_file);
    return unless $dump_fh;
    if ($committed) {
        $dump_fh->print("POST COMMIT UR CACHE INFO:\n");
    }
    elsif ($committed == 0) {
        $dump_fh->print("PRE COMMIT UR CACHE INFO:\n");
    }
    my %hash = ( 'Genome::TranscriptStructure' => scalar(keys %{$UR::Context::all_objects_loaded->{'Genome::TranscriptStructure'}}) );
    my @ss_keys = keys %{$UR::Context::all_objects_loaded->{'Genome::TranscriptStructure'}};

    my @ss_sample = map { $UR::Context::all_objects_loaded->{'Genome::TranscriptStructure'}->{$_}} @ss_keys[0..4];
    my $objects_loaded = $UR::Context::all_objects_cache_size;
    $dump_fh->print("all_objects_cache_size: $objects_loaded\n");
    $dump_fh->print(Dumper \%hash);
    $dump_fh->print("substructure samples:\n".Dumper \@ss_sample);
    $dump_fh->print("\n#########################################\n#######################################\n\n");
}

1;

# $Id$
