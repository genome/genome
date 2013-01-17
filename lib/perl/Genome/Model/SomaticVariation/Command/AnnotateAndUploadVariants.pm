package Genome::Model::SomaticVariation::Command::AnnotateAndUploadVariants;

use strict;
use warnings;
use Genome;
use Data::Dumper;

class Genome::Model::SomaticVariation::Command::AnnotateAndUploadVariants{
    is => 'Genome::Command::Base',
    has =>[
        build_id => {
            is => 'Integer',
            is_input => 1,
            is_output => 1,
            doc => 'build id of SomaticVariation model',
        },
        annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run   during the annotation step',
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3],
            is_optional => 1,
            is_input => 1,
        },
        build => {
            is => 'Genome::Model::Build::SomaticVariation',
            id_by => 'build_id',
        }
    ],
    has_param => [
        lsf_queue => {
            default => 'apipe',
        },
        lsf_resource => {
            default => "-M 20000000 -R 'select[mem>20000] rusage[mem=20000]'",
        },
        joinx_version => {
            is => "Text",
            doc => "Version of joinx to use",
            is_optional => 0,
            default => 1.6,
        },
        lsf_resource => { 
            is => 'Text',
            default => "-R 'rusage[mem=8000]' -M 8000000",
        },
    ],
};

sub execute{
    my $self = shift;
    my $build = $self->build;
    unless ($build){
        die $self->error_message("no build provided!");
    }

    $self->status_message("Executing Annotate and Upload step");

    my $version = 2;
    #my $version = GMT:BED:CONVERT::version();  TODO, something like this instead of hardcoding
    
    my %files;
    my %vcf_files;

    my ($tier1_snvs, $tier2_snvs, $tier1_indels, $tier2_indels);
    
    if ($build->snv_detection_strategy){
        $tier1_snvs = $build->data_set_path("effects/snvs.hq.novel.tier1",$version,'bed');
        unless(-e $tier1_snvs){
            die $self->error_message("No tier 1 snvs file for build!");
        }
        $files{'snvs.hq.tier1'} = $tier1_snvs;

        $tier2_snvs = $build->data_set_path("effects/snvs.hq.novel.tier2",$version,'bed');
        unless(-e $tier2_snvs){
            die $self->error_message("No tier 2 snvs file for build!");
        }
        $files{'snvs.hq.tier2'} = $tier2_snvs;
        my $snvs_vcf = $build->data_directory."/variants/snvs.vcf.gz";
        unless(-e $snvs_vcf) {
            die $self->error_message("No snvs vcf");
        }
        $vcf_files{'snvs.hq'} = $snvs_vcf;
    }

    if ($build->indel_detection_strategy){
        $tier1_indels = $build->data_set_path("effects/indels.hq.novel.tier1",$version,'bed');
        unless(-e $tier1_indels){
            die $self->error_message("No tier 1 indels file for build!");
        }
        $files{'indels.hq.tier1'} = $tier1_indels;

        $tier2_indels = $build->data_set_path("effects/indels.hq.novel.tier2",$version,'bed');
        unless(-e $tier2_indels){
            die $self->error_message("No tier 2 indels file for build!");
        }
        $files{'indels.hq.tier2'} = $tier2_indels;
    }

    #annotate variants
    my $annotator_version = $self->annotator_version;
    unless($annotator_version){
        die $self->error_message("No variant annotator version for build!");
    }

    my $annotation_build_id = $self->build->annotation_build->id;

    my %annotation_params = (
        use_version => $annotator_version,
        build_id => $annotation_build_id,
    );

    my $annotation_output_version=1;  #TODO, do we even have an annotation file format?  need to figure out how to resolve this w/ data_set_path
    my $upload_output_version = 1; #TODO, same issue here as with annotation output version


    for my $key (keys %files){
        my $variant_file = $files{$key};
        unless ($variant_file){
            $self->status_message("No detection strategy for $key. Skipping annotation and upload");
        }
        unless (-e $variant_file){
            die $self->error_message("File expected for annotating $key at $variant_file does not exist!  Failing");
        }

        my $annotated_file = $build->data_set_path("effects/$key",$annotation_output_version,'annotated');
        my $uploaded_file = $build->data_set_path("effects/$key", $upload_output_version, "uploaded");

        $annotation_params{variant_bed_file} = $variant_file;
        
        if (-s $variant_file){

            # Run annotation twice, once for the "none" annotation filter and once for the "top" annotation filter
            for my $annotation_filter ( qw(none top) ) {
                $annotation_params{annotation_filter} = $annotation_filter;

                # "None" will the be the default named output file, since annotated.none sounds a bit weird.
                my $final_output_file;
                if ($annotation_filter eq "none") {
                    $final_output_file = $annotated_file;
                    $annotation_params{output_file} = $final_output_file.".non-dedup";
                } 
                else {
                    $final_output_file = "$annotated_file.$annotation_filter";
                    $annotation_params{output_file} = "$final_output_file.non-dedup";
                }

                $self->status_message("Creating TranscriptVariants object with parameters: " . Dumper \%annotation_params);
                my $annotation_command = Genome::Model::Tools::Annotate::TranscriptVariants->create(%annotation_params);
                unless ($annotation_command){
                    die $self->error_message("Failed to create annotate command for $key. Params:\n".Data::Dumper::Dumper(\%annotation_params));
                }
                my $annotate_rv = $annotation_command->execute;
                my $annotate_err = $@;
                unless($annotate_rv){
                    die $self->error_message("Failed to execute annotate command for $key(err: $annotate_err) from params:\n" . Data::Dumper::Dumper(\%annotation_params));
                }
                unless(-s $annotation_params{output_file}){
                    die $self->error_message("No output from annotate command for $key. Params:\n" . Data::Dumper::Dumper(\%annotation_params));
                }
                #Deduplicate the annotation file
                my $in = Genome::Sys->open_file_for_reading($annotation_params{output_file});
                my $output = Genome::Sys->open_file_for_writing($final_output_file);
                my %seen;
                while(my $line = <$in>) {
                    chomp $line;
                    if ($seen{$line}) {
                        next;
                    }
                    $seen{$line} = 1;
                    print $output $line."\n";
                }
                $in->close;
                $output->close;
                my $rm_cmd = "rm -f ".$annotation_params{output_file};
                `$rm_cmd`;
                
                #Make a version of the annotation file without a header
                `mv $final_output_file $final_output_file.header`;
                `grep -v "^chromosome_name" $final_output_file.header > $final_output_file`;
            }

            #upload variants
            my %upload_params = (
                build_id => $self->build_id, 
            );

            $upload_params{variant_file} = $variant_file;
            $upload_params{annotation_file} = $annotated_file;
            $upload_params{output_file} = $uploaded_file;

            my $upload_command = Genome::Model::Tools::Somatic::UploadVariants->create(%upload_params);
            unless ($upload_command){
                die $self->error_message("Failed to create upload command for $key. Params:\n". Data::Dumper::Dumper(\%upload_params));
            }
            #my $upload_rv = $upload_command->execute;
            my $upload_rv =  1;  #TODO turn this on
            my $upload_err = $@;
            unless ($upload_rv){
                die $self->error_message("Failed to execute upload command for $key (err: $upload_err). Params:\n".Dumper(\%upload_params));
            }
            unless(-s $upload_params{output_file} or 1){
                die $self->error_message("No output from upload command for $key. Params:\n" . Data::Dumper::Dumper(\%upload_params));
            }
        }
        else{
            $self->status_message("No variants present for $key, skipping annotation and upload");
            File::Copy::copy($variant_file, $annotated_file);
            File::Copy::copy($variant_file, "$annotated_file.top");
            File::Copy::copy($variant_file, $uploaded_file);
        }
    }

    #Annotate vcfs with dbsnp ids
    $self->status_message("Adding dbsnp ids to vcf");
    
    if ($build->previously_discovered_variations_build) {
        my $annotation_vcf = $build->previously_discovered_variations_build->snvs_vcf;
        if (-e $annotation_vcf) {
            for my $key (keys(%vcf_files)) {
                my $variant_file = $vcf_files{$key};
                my $output_file = $variant_file.".annotated.vcf.gz";
                $output_file =~ s/vcf.gz.//;
                my $info_string = $build->processing_profile->vcf_annotate_dbsnp_info_field_string eq "NO_INFO" ? "" : $build->processing_profile->vcf_annotate_dbsnp_info_field_string;
                my $info = $info_string eq "" ? 0 : 1;
                my $vcf_annotator = Genome::Model::Tools::Joinx::VcfAnnotate->execute(
                    input_file=> $variant_file,
                    annotation_file=>$annotation_vcf,
                    output_file=>$output_file,
                    use_bgzip=>1,
                    info_fields=>$info_string,
                    info => $info,
                    use_version => $self->joinx_version,
                ) || die "Failed to execute Joinx Vcf annotation using db: $annotation_vcf";
                $self->status_message("Successfully annotated VCF with information from $annotation_vcf");
            }
            foreach my $annotation_prefix ("snvs.hq.tier1", "snvs.hq.tier2") {
                my $annotation_top_file = $build->data_set_path("effects/$annotation_prefix", $annotation_output_version, "annotated.top");
                next unless (-e $annotation_top_file);
                my $vcf_file = $vcf_files{"snvs.hq"};
                $vcf_file = $vcf_file.".annotated.vcf.gz";
                $vcf_file =~ s/vcf.gz.//;
                $self->status_message("Adding rsid column for $annotation_prefix from vcf file ".$vcf_file);
                my $rv = Genome::Model::Tools::Annotate::AddRsid->execute(
                    anno_file => $annotation_top_file,
                    vcf_file => $vcf_file,
                    output_file => $build->data_set_path("effects/$annotation_prefix.rsid", $annotation_output_version, "annotated.top"),
                );
            }
        }
        else {
            $self->warning_message("No snvs vcf available for previously_discovered_variations_build");
        }
    }
    
    #Annotate SVs. Need to handle mouse vs human and human 36 vs 37
    $self->status_message("Executing sv annotation");

    if ($build->sv_detection_strategy) {
        my $sv_sr  = $build->final_result_for_variant_type('sv');
        if ($sv_sr) {
            my $sv_dir = $sv_sr->output_dir;
            if ($sv_dir and -d $sv_dir) {
                my $sv_file = $sv_dir . '/svs.merge.file.somatic';
                if ($sv_file and -s $sv_file) {
                    my $species_name = $build->subject->species_name;
                    if ($species_name and $species_name =~/^(human|mouse)$/) {
                        my $sv_annot_out_file = $build->data_directory . '/effects/svs.hq.annotated';
                        my $fusion_out_file   = $build->data_directory . '/effects/svs.hq.fusion_transcripts.out';
                        my @annotator_list    = qw(Transcripts FusionTranscripts);

                        my $cancer_gene_list = join("/",Genome::Sys->dbpath(join("/","cancer-gene-list",$species_name), 1),"Cancer_genes.csv");

                        my %params = (
                            input_file  => $sv_file,
                            output_file => $sv_annot_out_file,
                            fusion_transcripts_fusion_output_file => $fusion_out_file,
                            annotation_build_id => $annotation_build_id,
                            annotator_list      => \@annotator_list,
                            transcripts_print_flanking_genes => 1,
                            transcripts_cancer_gene_list     => $cancer_gene_list,
                            chrA_column       => 1,
                            bpA_column        => 2,
                            chrB_column       => 4,
                            bpB_column        => 5,
                            event_type_column => 7,
                            orient_column     => 8,
                            score_column      => 12,  
                        );
                        
                        my $ref_seq_build = $build->reference_sequence_build;
                        if ($species_name eq 'human') {
                            my $annot_file = $self->_get_human_sv_annot_file($ref_seq_build);

                            if ($annot_file) {
                                push @annotator_list, 'Dbsnp', 'Segdup', 'Dbvar';
                                %params = (
                                    %params,
                                    annotator_list => \@annotator_list,
                                    dbsnp_annotation_file  => $annot_file->{dbsnp},
                                    segdup_annotation_file => $annot_file->{segdup},
                                    dbvar_annotation_file  => $annot_file->{dbvar},
                                );
                            }
                            else {
                                $self->warning_message('Human annotation for dbsnp, dbvar and segdup is only set for build 36 and 37');
                            }
                        }
                        elsif ($species_name eq 'mouse') {
                            unless ($ref_seq_build->id == 107494762) {#UCSC-mouse build37, mouse ref seq used in production
                                $self->warning_message('For now sv annotation only run on UCSC-mouse build37');
                                %params = ();
                            }
                        }

                        if (%params) {
                            my $sv_annot = Genome::Model::Tools::Annotate::Sv::Combine->create(%params);
                            my $rv = $sv_annot->execute;    
                            $self->warning_message("sv annotation did not finish ok") unless $rv == 1;
                        }
                    }
                    else {
                        $self->warning_message('The build is not for human or mouse. Skip sv annotation');
                    }
                }
                else {
                    $self->warning_message("combine sv somatic file: $sv_file is not valid");
                }
            }
            else {
                $self->warning_message("sv dir: $sv_dir is not valid");
            }
        }
        else {
            $self->warning_message('No sv software result for this build');
        }
    }

    my $annovar = Genome::Model::Tools::Annovar::AnnotateVariation->execute(
        input_file => $build->data_directory."/variants/snvs.annotated.vcf.gz",
        buildver => "hg19",
        table_names => ["wgEncodeRegDnaseClustered", "wgEncodeRegTfbsClustered", "bed"],
        outfile => $build->data_directory."/effects/snvs.annovar",
        annotation_type => "regionanno",
        scorecolumn => 5,
        bedfiles => ["/gscuser/aregier/scratch/dhs_promoters/gene_names2.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/segway.hg19.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/wgEncodeBroadHmmGm12878HMM.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/wgEncodeBroadHmmH1hescHMM.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/wgEncodeBroadHmmHepg2HMM.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/wgEncodeBroadHmmHmecHMM.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/wgEncodeBroadHmmHsmmHMM.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/wgEncodeBroadHmmHuvecHMM.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/wgEncodeBroadHmmK562HMM.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/wgEncodeBroadHmmNhekHMM.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/wgEncodeBroadHmmNhlfHMM.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/DRM_transcript_pairs.bed",
                     "/gscuser/aregier/scratch/dhs_promoters/Open_Chromatin_Supp_Table_4.bed"],
    );
    unless ($annovar) {
        $self->error_message("Annovar failed");
        return;
    }
    if ($build->data_set_path("effects/snvs.hq.tier2", $annotation_output_version, "annotated.top.header")){
        my $temp = Genome::Sys->create_temp_file_path;
        my $counter = 0;
        my $in_file = $build->data_set_path("effects/snvs.hq.tier2", $annotation_output_version, "annotated.top.header");
        foreach my $table_name (("wgEncodeRegDnaseClustered", "wgEncodeRegTfbsClustered", "bed_gene_names2", "bed_segway.hg19", "bed_wgEncodeBroadHmmGm12878HMM", "bed_wgEncodeBroadHmmH1hescHMM", "bed_wgEncodeBroadHmmHepg2HMM", "bed_wgEncodeBroadHmmHmecHMM", "bed_wgEncodeBroadHmmHsmmHMM", "bed_wgEncodeBroadHmmHuvecHMM", "bed_wgEncodeBroadHmmK562HMM", "bed_wgEncodeBroadHmmNhekHMM", "bed_wgEncodeBroadHmmNhlfHMM", "bed_DRM_transcript_pairs", "bed_Open_Chromatin_Supp_Table_4")) {
            $counter++;
            my $append = Genome::Model::Tools::Annotate::AppendColumns->execute(
                additional_columns_file => $build->data_directory."/effects/snvs.annovar.hg19_".$table_name,
                input_variants => $in_file,
                output_file => "$temp.$counter",
                column_to_append => 2,
                header => $table_name,
                chrom_column => 3,
                start_column => 4,
                stop_column => 5,
            );
            unless ($append) {
                $self->error_message("Append columns failed for tier2 snvs");
                return;
            }
            $in_file = "$temp.$counter";
        }
        my $cmd = "mv $temp.$counter ".$build->data_set_path("effects/snvs.hq.tier2", $annotation_output_version, "annotated.top.annovar");
        `$cmd`;
    }
    if ($build->data_set_path("effects/snvs.hq.novel.tier3", $version, "bed")) {
        my $temp_file = Genome::Sys->create_temp_file_path;
        my $convert = Genome::Model::Tools::Bed::Convert::BedToAnnotation->execute(
            snv_file => $build->data_set_path("effects/snvs.hq.novel.tier3", $version, "bed"),
            output => $temp_file,
            annotator_version => 3,
        );
        #Add a header
        my $output_file = Genome::Sys->open_file_for_writing($build->data_set_path("effects/snvs.hq.novel.tier3", $annotation_output_version, "converted-anno"));
        $output_file->print("chromosome_name\tstart\tstop\treference\tvariant\n");
        my $in = Genome::Sys->open_file_for_reading($temp_file);
        while(my $line = <$in>) {
            $output_file->print($line);
        }
        $in->close;
        $output_file->close;

        unless ($convert) {
            $self->error_message("Conversion from bed to anno coords failed for tier3 snvs");
            return;
        }
        my $temp = Genome::Sys->create_temp_file_path;
        my $in_file = $build->data_set_path("effects/snvs.hq.novel.tier3", $annotation_output_version, "converted-anno");
        my $counter = 0;
        foreach my $table_name (("wgEncodeRegDnaseClustered", "wgEncodeRegTfbsClustered", "bed_gene_names2", "bed_segway.hg19", "bed_wgEncodeBroadHmmGm12878HMM", "bed_wgEncodeBroadHmmH1hescHMM", "bed_wgEncodeBroadHmmHepg2HMM", "bed_wgEncodeBroadHmmHmecHMM", "bed_wgEncodeBroadHmmHsmmHMM", "bed_wgEncodeBroadHmmHuvecHMM", "bed_wgEncodeBroadHmmK562HMM", "bed_wgEncodeBroadHmmNhekHMM", "bed_wgEncodeBroadHmmNhlfHMM", "bed_DRM_transcript_pairs", "bed_Open_Chromatin_Supp_Table_4")) {
            $counter++;
            my $append = Genome::Model::Tools::Annotate::AppendColumns->execute(
                additional_columns_file => $build->data_directory."/effects/snvs.annovar.hg19_".$table_name,
                input_variants => $in_file,
                output_file => "$temp.$counter",
                column_to_append => 2,
                header => $table_name,
                chrom_column => 3,
                start_column => 4,
                stop_column => 5,
            );
            unless($append) {
                $self->error_message("Append columns failed for tier3 snvs");
                return;
            }
            $in_file = "$temp.$counter";
        }
        my $cmd = "mv $temp.$counter ".$build->data_set_path("effects/snvs.hq.novel.tier3", $annotation_output_version, "converted-anno.annovar");
        `$cmd`;
    }

    #upload variants

    ##upload metrics
    $self->status_message("Uploading metrics");
    eval{
        my ($tier3_snvs, $tier4_snvs, $tier3_indels, $tier4_indels, $svs);
        $tier3_snvs = $build->data_set_path("effects/snvs.hq.novel.tier3",$version,'bed');
        $tier4_snvs = $build->data_set_path("effects/snvs.hq.novel.tier4",$version,'bed');
        $tier3_indels = $build->data_set_path("effects/indels.hq.novel.tier3",$version,'bed');
        $tier4_indels = $build->data_set_path("effects/indels.hq.novel.tier4",$version,'bed');

        $svs = $build->data_set_path("variants/svs.hq",$version,'bed');

        my ($number_of_tier1_snvs, $number_of_tier2_snvs, $number_of_tier3_snvs, $number_of_tier4_snvs, $number_of_tier1_indels, $number_of_tier2_indels, $number_of_tier3_indels, $number_of_tier4_indels, $number_of_svs);
        if($build->snv_detection_strategy) {
            chomp($number_of_tier1_snvs = `wc -l $tier1_snvs | cut -f 1 -d' '`);
            chomp($number_of_tier2_snvs = `wc -l $tier2_snvs | cut -f 1 -d' '`);
            chomp($number_of_tier3_snvs = `wc -l $tier3_snvs | cut -f 1 -d' '`);
            chomp($number_of_tier4_snvs = `wc -l $tier4_snvs | cut -f 1 -d' '`);
            $build->set_metric("tier1_snv_count", $number_of_tier1_snvs);
            $build->set_metric("tier2_snv_count", $number_of_tier2_snvs);
            $build->set_metric("tier3_snv_count", $number_of_tier3_snvs);
            $build->set_metric("tier4_snv_count", $number_of_tier4_snvs);
        }


        if($build->indel_detection_strategy) {
            chomp($number_of_tier1_indels = `wc -l $tier1_indels | cut -f 1 -d' '`);
            chomp($number_of_tier2_indels = `wc -l $tier2_indels | cut -f 1 -d' '`);
            chomp($number_of_tier3_indels = `wc -l $tier3_indels | cut -f 1 -d' '`);
            chomp($number_of_tier4_indels = `wc -l $tier4_indels | cut -f 1 -d' '`);
            $build->set_metric("tier1_indel_count", $number_of_tier1_indels);
            $build->set_metric("tier2_indel_count", $number_of_tier2_indels);
            $build->set_metric("tier3_indel_count", $number_of_tier3_indels);
            $build->set_metric("tier4_indel_count", $number_of_tier4_indels);
        }

        if($build->sv_detection_strategy) {
            chomp($number_of_svs = `wc -l $svs | cut -f 1 -d' '`); 
            $build->set_metric("sv_count", $number_of_svs);
        }
    }; 
    $self->status_message("Metric setting problem: $@") if $@;

    return 1;
}

#only run on human build 36 and 37 for now
sub _get_human_sv_annot_file {
    my ($self, $ref_seq_build) = @_;
    my $version = $ref_seq_build->version;
    my $type;
    my $dbsnp_version;
    my $species = "human";

    if ($version) {
        if ($version =~ /36/) {
            $dbsnp_version = "130";
            $type = "build36";
        }
        elsif ($version =~ /37/) {
            $dbsnp_version = "132";
            $type = "build37";
        }
    }
    
    if ($ref_seq_build->id == 109104543) {
        $type = "build36";
        $dbsnp_version = "130";
    }
    return unless $type and $dbsnp_version;
    $self->status_message("Getting annotation files for species $species of version $type, dbsnp version $dbsnp_version");
    my $segdup = Genome::Db->get(source_name => 'ucsc', database_name => $species, external_version => $type);
    my $dbsnp = Genome::Db->get(source_name => 'genome-db-dbsnp', database_name => "$species/$type", external_version => $dbsnp_version);
    my $dbvar = Genome::Db->get(source_name => 'dbvar', database_name => $species, external_version => $type);
    return {
        segdup => $segdup->data_directory . '/segdup.tsv',
        dbsnp  => $dbsnp->data_directory . '/dbsnp.csv',
        dbvar  => $dbvar->data_directory . '/dbvar.tsv',
        };
}

1;
