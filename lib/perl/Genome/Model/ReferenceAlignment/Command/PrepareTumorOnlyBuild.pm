package Genome::Model::ReferenceAlignment::Command::PrepareTumorOnlyBuild;

class Genome::Model::ReferenceAlignment::Command::PrepareTumorOnlyBuild {
    is => 'Command::V2',
    doc => 'Add dbsnp filtering for results of refalign build on tumor only model',
    has => [
        build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            shell_args_position => 1,
            doc => 'build to prepare',
        },
        dbsnp_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'dbSNP build to filter on, if not given, will use dbSNP build from the build input',
            is_optional => 1,
        },
        bed_version => {
            is => 'Text',
            doc => 'Version of bed file to require.',
            default => '2',
        },
        previously_detected_snvs => {
            is => 'Text',
            doc => 'Intersect snvs with this list to filter out previously detected snvs',
            is_optional => 1,
        },
        move_existing => {
            is => 'Boolean',
            doc => 'Set this to move existing ./effects directories out of the way before proceding.',
            default => 0,
        },
    ],
};

sub help_detail {
    return <<EOS
PrepareTumorOnlyBuild is a script which is meant to be run on a reference-alignment build with a completed
detect-variants step. This script will pull the snvs and indels then create an effects subdir on the build's
data-directory. This script then splits the snvs into dbSNP and novel, based on either the dbsnp_build param
or the dbsnp_build param on the refalign build passed in.
EOS
}
#'

sub execute {
    my $self = shift;

    #reference alignment build object
    my $build = $self->build;

    #get the appropiate dbSNP build
    my $dbsnp_build = defined($self->dbsnp_build) ? $self->dbsnp_build : $build->dbsnp_build;
    my $dbsnp_file =  $dbsnp_build->snvs_bed;
    unless(defined($dbsnp_build)){
        $self->status_message("No dbsnp build found, nothing to do. Exiting.");
        return 1;
    }

    my $build_dir = $build->data_directory;
    my $version = $self->bed_version;
    my $snv_bed_file = $build_dir."/variants/snvs.hq.v".$version.".bed";

    unless(-e $snv_bed_file){
        die $self->error_message("Could not locate snvs file at: ".$snv_bed_file."\n");
    }

    my $output_dir = $build_dir."/effects";
    my $in_roi = $output_dir."/snvs.hq.in_roi.v".$version.".bed";
    my $not_in_roi = $output_dir."/snvs.hq.not_in_roi.v".$version.".bed";
    my $in_dbsnp_file = $output_dir."/snvs.hq.in_dbsnp.v".$version.".bed";
    my $not_in_dbsnp_file = $output_dir."/snvs.hq.not_in_dbsnp.v".$version.".bed";
    my $previously_detected_file = $output_dir."/snvs.hq.previously_detected.v".$version.".bed";
    my $novel_file = $output_dir."/snvs.hq.novel.v".$version.".bed";

    my $annotation_build = $build->annotation_reference_build;
    my $tier_file_location = $annotation_build->tiering_bed_files_by_version(3);
    unless(defined($tier_file_location)){
        die $self->error_message("Could not locate tiering files!");
    }

    #create the output directory
    unless(-d $output_dir){
        Genome::Sys->create_directory($output_dir);
    }
    elsif ($self->move_existing) {
        File::Copy::move($output_dir, $output_dir."_previous");
        Genome::Sys->create_directory($output_dir);
    }

    #my $roi_file = $output_dir."/target_region_bid_".$build->id.".bed";
    unless(-e $in_roi && -e $not_in_roi){
        my $roi_name = $build->model->region_of_interest_set_name;
        my $roi = Genome::FeatureList->get( name => $roi_name);
        my @roi = split /\n/, $roi->processed_bed_file_content;
        my $roi_temp = Genome::Sys->create_temp_file_path;
        my $rtfh = Genome::Sys->open_file_for_writing($roi_temp);
        for my $line (@roi){
            my @line = split /\s/, $line;
            print $rtfh join("\t",($line[0],$line[1],$line[2]))."\n";
        }
        $rtfh->close;

        my $roi_intersection_cmd = Genome::Model::Tools::Joinx::Intersect->create(
            input_file_a => $snv_bed_file,
            input_file_b => $roi_temp,
            miss_a_file => $not_in_roi,
            output_file => $in_roi,
        );
        unless($roi_intersection_cmd->execute){
            die $self->error_message("Could not complete ROI intersection.");
        }
        $self->status_message("Completed ROI intersection: ".$not_in_roi."\n");
    }
    #run dbSNP intersection
    unless(-e $not_in_dbsnp_file && -e $in_dbsnp_file){
        my $dbsnp_intersection = Genome::Model::Tools::Joinx::Intersect->create(
            input_file_a => $in_roi,
            input_file_b => $dbsnp_file,
            miss_a_file => $not_in_dbsnp_file,
            output_file => $in_dbsnp_file,
            dbsnp_match => 1,
        );
        unless ($dbsnp_intersection){
            die $self->error_message("Couldn't create joinx intersection tool!");
        }
        my $snv_rv = $dbsnp_intersection->execute();
        my $snv_err = $@;
        unless ($snv_rv){
            die $self->error_message("Failed to execute joinx intersection (err: $snv_err )");
        }
    }

    my @snv_tiers;
    my $snv_files_present=1;

    #intersect with previously detected snvs
    my $previously_detected_snvs = $self->previously_detected_snvs;
    if(defined($previously_detected_snvs)){
        unless(-e $novel_file && -e $previously_detected_file){
            my $dbsnp_intersection = Genome::Model::Tools::Joinx::Intersect->create(
                input_file_a => $not_in_dbsnp_file,
                input_file_b => $previously_detected_snvs,
                miss_a_file => $novel_file,
                output_file => $previously_detected_file,
                );
            unless ($dbsnp_intersection){
                die $self->error_message("Couldn't create joinx intersection tool!");
            }
            my $snv_rv = $dbsnp_intersection->execute();
            my $snv_err = $@;
            unless ($snv_rv){
                die $self->error_message("Failed to execute joinx intersection (err: $snv_err )");
            }
        }

        #compose snv and indel tiering output file names
        @snv_tiers = map{ "snvs.hq.novel.tier".$_.".v".$version.".bed"} (1..4);
        for (1..4) {
            my $f = -e $output_dir."/".$snv_tiers[$_];
            $snv_files_present = $snv_files_present && $f;
        }

    } else { #no prev file, tier and annotated the not_in_dbsnp file

        #compose snv and indel tiering output file names
        @snv_tiers = map{ "snvs.hq.not_in_dbsnp.tier".$_.".v".$version.".bed"} (1..4);
        for (1..4) {
            my $f = -e $output_dir."/".$snv_tiers[$_];
            $snv_files_present = $snv_files_present && $f;
        }
        $novel_file = $not_in_dbsnp_file
    }



    #tier snvs
    unless( $snv_files_present ){
        my $snv_tier_cmd = Genome::Model::Tools::FastTier::FastTier->create(
            variant_bed_file => $novel_file,
            tier_file_location => $tier_file_location,
            tier1_output => $output_dir."/".$snv_tiers[0],
            tier2_output => $output_dir."/".$snv_tiers[1],
            tier3_output => $output_dir."/".$snv_tiers[2],
            tier4_output => $output_dir."/".$snv_tiers[3],
        );
        unless($snv_tier_cmd){
            die $self->error_message("Could not create tiering command for snvs!");
        }
        unless($snv_tier_cmd->execute){
            die $self->error_message("Tiering Snvs did not succeed!");
        }
    }

    #prepare annotation output file names
    my $snv_tier1_anno_output = $output_dir."/snvs.hq.novel.tier1.annotated";
    my $snv_tier2_anno_output = $output_dir."/snvs.hq.novel.tier2.annotated";
    if(!(defined($previously_detected_snvs))){
        $snv_tier1_anno_output = "snvs.hq.not_in_dbsnp.tier1.annotated";
        $snv_tier1_anno_output = "snvs.hq.not_in_dbsnp.tier2.annotated";
    }

    #Annotate Snvs
    unless(-e $snv_tier1_anno_output && -e $snv_tier2_anno_output){
        my $snv_tier1_anno_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
            variant_bed_file => $output_dir."/".$snv_tiers[0],
            output_file => $snv_tier1_anno_output,
            annotation_filter => "top",
            build_id => $annotation_build->id,
            use_version => 2,
        );
        unless($snv_tier1_anno_cmd){
            die $self->error_message("Could not create snv_tier1_anno_cmd!");
        }
        unless($snv_tier1_anno_cmd->execute){
            die $self->error_message("Could not complete snv_tier1_anno_cmd");
        }
        my $snv_tier2_anno_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
            variant_bed_file => $output_dir."/".$snv_tiers[1],
            output_file => $snv_tier2_anno_output,
            annotation_filter => "top",
            build_id => $annotation_build->id,
            use_version => 2,
        );
        unless($snv_tier2_anno_cmd){
            die $self->error_message("Could not create snv_tier1_anno_cmd!");
        }
        unless($snv_tier2_anno_cmd->execute){
            die $self->error_message("Could not complete snv_tier1_anno_cmd");
        }
    }

    #reallocate build's data_directory since we have added some data
    my $build_allocation = $build->disk_allocation;
    $build_allocation->reallocate;
    return 1;
}

1;
