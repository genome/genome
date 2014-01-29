package Genome::Model::SomaticValidation::Command::ReviewVariants;

use warnings;
use strict;
use IO::File;
use Genome;
use Sort::Naturally qw(nsort);
use Genome::Info::IUB;
use File::Basename;

class Genome::Model::SomaticValidation::Command::ReviewVariants {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            doc => 'the SomaticValidation build whose results to compare to the candidate variants',
            id_by => 'build_id',
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            via => 'build',
            to => 'reference_sequence_build',
            is_optional => 1, #FIXME not really, but this makes the CLI happy
        },
        somatic_variation_model => {
            is => 'Genome::Model::SomaticVariation',
            is_optional => 1,
            id_by => 'somatic_variation_model_id',
        },
    ],

    has_input => [
        build_id => {
            is => 'Text',
            doc => "ID of SomaticValidation build",
            is_output => 1,
        },
        output_dir => {
            is => 'Text',
            is_optional => 1,
            doc => "Directory where output will be stored (under a subdirectory with the sample name) [defaults to BUILD_DIRECTORY/validation]",
        },
    ],

    has_optional_input => [
        process_indels => {
            is => 'Boolean',
            doc => 'Indicate whether indels will be processed for review',
            default_value => 0,
        },

        variant_list =>{
            is => 'Text',
            is_optional => 1,
            doc => "list of variants that we're trying to validate, in bed format. If not specified, it will try to grab this list from the model. If no list is found, it will assume that this is an extension experiment (and that there are no variants to validate)",
        },

        somatic_variation_model_id =>{
            is => 'Text',
            is_optional => 1,
            doc => "somatic variation model that was used to call the variants in the first place. Only used if --review-non-validated-calls is true. If not specified, we'll try to grab this from the som-val model link",
        },

        igv_reference_name =>{
            is => 'Text',
            is_optional => 1,
            doc => "name of the igv reference to use",
            is_calculated => 1,
            calculate_from => ['reference_sequence_build'],
            calculate => q{
                if($reference_sequence_build->species_name eq "human") {
                    if($reference_sequence_build->name eq 'NCBI-human-build36') {
                        return 'reference_build36';
                    } else {
                        return 'b' . $reference_sequence_build->version;
                    }
                } else {
                    return;
                }
            },
            example_values => ["reference_build36"],
        },

        filter_sites =>{
            is => 'Text',
            is_optional => 1,
            doc => "comma separated list of bed files containing sites to be removed (example - removing cell-line sites from tumors grown on them). Only sites with these exact coordinates and that match the ref and var alleles are removed.",
        },

        filter_regions =>{
            is => 'Text',
            is_optional => 1,
            doc => "comma separated list of bed files containing regions to be removed (example - removing paralog regions) Any sites intersecting these regions are removed.",
        },

        restrict_to_target_regions =>{
            is => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => "only keep snv calls within the target regions. These are pulled from the build",
        },

        tier1_only =>{
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => "only keep and review calls that are tier 1",
        },

        tier_file_location =>{
            is => 'String',
            is_optional => 1,
            is_calculated => 1,
            calculate_from => ['build'],
            calculate => q { $build->annotation_build->tiering_bed_files_by_version($build->processing_profile->tiering_version) },
            doc => "if tier1-only is specified, this needs to be a path to the appropriate tiering files",
        },

        recover_non_validated_calls =>{
            is => 'Boolean',
            is_optional => 1,
            doc => "Do readcounting on calls that are non-validated, find sites that have high coverage (>20x tumor and normal) in the original, but poor (<8x tumor/ <6x norm) coverage in the validation. Send these for review",
            default => 0,
        },

        variant_list_is_one_based => {
            is => 'Boolean',
            is_optional => 1,
            doc => "The variant list you provide is in annotation (one-based) format, instead of bed. This flag fixes that.",
            default => 0,
        },

        dbsnp_filter => {
            is => 'String',
            is_optional => 1,
            doc => "path to a dbsnp bed file of sites that should be removed from consideration",
        },

        ##restrict to targeted region - grab from build

        # read_review => {
        #     is => 'Boolean',
        #     doc => "Read existing manual review files and create WU annotation files per case",
        #     is_optional => 1,
        #     default => 0
        #   },

    ],
};

sub sub_command_category { 'pipeline steps' }

sub help_detail {
    return <<HELP;
Given a SomaticValidation model, this tool will gather the resulting variants, remove
off-target sites, tier the variants, optionally filter them, and match them up with the
initial predictions sent for validation.  It will then divide them into categories (validated,
non-validated, and new calls). New calls are prepped for manual review in the review/ directory.
HELP
}

sub _doc_authors {
    return <<AUTHS;
 Chris Miller
AUTHS
}

sub bedToAnno{
    my ($chr,$start,$stop,$ref,$var) = split("\t",$_[0]);
    #print STDERR join("|",($chr,$start,$stop,$ref,$var)) . "\n";
    if ($ref =~ /^\-/){ #indel INS
        $stop = $stop+1;
    } else { #indel DEL or SNV
        $start = $start+1;
    }
    return(join("\t",($chr,$start,$stop,$ref,$var)));
}

sub annoToBed{
    my ($chr,$start,$stop,$ref,$var) = split("\t",$_[0]);
    if ($ref =~ /^\-/){ #indel INS
        $stop = $stop-1;
    } else { #indel DEL or SNV
        $start = $start-1;
    }
    return(join("\t",($chr,$start,$stop,$ref,$var)));
}

sub intersects{
    my ($st,$sp,$st2,$sp2) = @_;
    if((($sp2 >= $st) && ($sp2 <= $sp)) ||
        (($sp >= $st2) && ($sp <= $sp2))){
        return 1;
    }
    return 0;
}


sub execute {
    my $self = shift;

    my $build = $self->build;
    my $tumor_only = not($build->normal_sample);
    my $output_dir = $self->output_dir;

    if($output_dir) {
        $output_dir =~ s/(\/)+$//; # Remove trailing forward-slashes if any
        unless(-e $output_dir) {
            die $self->error_message('The specified output directory does not exist.');
        }
    } else {
        $output_dir = join('/', $build->data_directory, 'validation');
        Genome::Sys->create_directory($output_dir);
    }
    $self->output_dir($output_dir);

    #FIXME This should probably grab specific tumor_build and normal_build of the specified build.
    unless($self->somatic_variation_model){
        if(defined($build->snv_variant_list)){
            if(defined($build->snv_variant_list->source_build)){
                $self->somatic_variation_model($build->snv_variant_list->source_build->model);
            }
        }
    }

    #grab the info from all of the models
    my %bams; # Hash to store the model info

    #my $ref_seq_build = $build->reference_sequence_build;
    #my $ref_seq_fasta = $ref_seq_build->full_consensus_path('fa');
    my $sample_name = $build->tumor_sample->name;
    $self->debug_message('processing build ' . $build->id . ' with sample_name: ' . $sample_name);
    #my $tumor_bam = $build->tumor_bam;
    my $build_dir = $build->data_directory;

    #my $normal_bam;
    #unless($tumor_only){
    #    $normal_bam = $build->normal_bam;
    #}

    unless ($build->snv_detection_strategy) {
        $self->warning_message("snv_detection_strategy undefined, skipping Review Variants. This tool could be rewritten to allow for no snv_detection_strategy if a variant list is provided.");
        return 1;
    }

    # Ensure the necessary files exist in this build
    my $snv_file = "$build_dir/variants/snvs.hq.bed";
    unless( -e $snv_file ){
        die $self->error_message("SNV file not found at $snv_file");
    }


    my $indel_file = "$build_dir/variants/indels.hq.bed";
    unless(not $self->process_indels or -e $indel_file ){
        die $self->error_message("INDEL file not found at $indel_file");
    }

    # create subdirectories, get files in place
    for my $subdir ("snvs", ($self->process_indels? "indels" : ()), "review", "snvs/readcounts") {
        Genome::Sys->create_directory(join("/", $output_dir, $subdir));
    }

    #generate variant files with any duplicates removed
    my $clean_snv_file = $self->generate_clean_file($snv_file);
    my $clean_indel_file = $self->generate_clean_file($indel_file) if $self->process_indels;

    #store the previously called variants into a hash
    my $prev_calls = $self->parse_prior_variants;

    for my $f ($clean_snv_file, ($self->process_indels? $clean_indel_file : ())) {
        $self->classify_variants($f, $prev_calls);
    }

    my $new_indel_file = "$output_dir/indels/indels.newcalls";
    my $new_snv_file = "$output_dir/snvs/snvs.newcalls";

    #from here on out there are several optional processes.  Keep a reference to the most completely processed file in this variable
    my $current_snv_file = $new_snv_file;
    my $current_indel_file = $new_indel_file;

    #-------------------------------------------------
    #filter out the off-target regions, if target regions are available

    if($self->restrict_to_target_regions and $build->region_of_interest_set){
        my $target_regions = $self->prepare_roi($build->region_of_interest_set);

        $current_snv_file = $self->region_limit_variants($new_snv_file, $target_regions);
        $current_indel_file = $self->region_limit_variants($new_indel_file, $target_regions) if $self->process_indels;
    }

    ##------------------------------------------------------
    #remove all but tier 1 sites, if that option is specified
    if($self->tier1_only){
        $self->debug_message("Doing Tiering...");
        $current_snv_file = $self->tier_variant_file($current_snv_file);
        $current_indel_file = $self->tier_variant_file($current_indel_file) if $self->process_indels;
    }

    ##-------------------------------------------------
    #use joinx to remove dbsnp sites
    if(defined($self->dbsnp_filter)){
        $self->debug_message("Applying dbsnp filter...");

        $current_snv_file = $self->dbsnp_filter_variant_file($current_snv_file);
        $current_indel_file = $self->dbsnp_filter_variant_file($current_indel_file) if $self->process_indels;
    }

    #-------------------------------------------------
    #remove filter sites specified by the user
    if(defined($self->filter_sites)){
        $self->debug_message('Applying user-supplied filter...');

        my $filter_sites = $self->prepare_user_filter_sites($self->filter_sites);

        $current_snv_file = $self->filter_variants($filter_sites, $current_snv_file);
        $current_indel_file = $self->filter_variants($filter_sites, $current_indel_file) if $self->process_indels;
    }

    #-------------------------------------------------
    #remove filter regions specified by the user
    if(defined($self->filter_regions)){
        my $filter_region_file = $self->prepare_user_filter_regions($self->filter_regions);

        $self->debug_message("Removing user-specified filter...");

        $current_snv_file = $self->_filter_regions($filter_region_file, $current_snv_file);
        $current_indel_file = $self->_filter_regions($filter_region_file, $current_indel_file) if $self->process_indels;
    }


    my $notvalidated_snvs_file = "$output_dir/snvs/snvs.notvalidated";
    my $validated_snvs_file = "$output_dir/snvs/snvs.validated";
    #add readcounts
    $self->debug_message("Getting readcounts...");
    foreach my $file ($validated_snvs_file,$notvalidated_snvs_file,$current_snv_file){
        if( -s $file ){
            $self->generate_readcounts($file);
        }
    }



    ## we're not going to do this in most cases - if not called in validation, we just ditch the call,
    ## but it may be useful in some cases.

    #-------------------------------------------------
    #look at the calls that were missed (case 3 - called in original, failed validation)
    #to determine whether they were missed due to poor coverage.
    #if coverage is fine, dump them (most), but if coverage was poor in validation (and good in wgs), send for review
    #we can only really do this for snvs at the moment, until indel readcounting is tweaked
    if(-s $notvalidated_snvs_file and $self->recover_non_validated_calls) {
        $self->recover_low_coverage_variants($notvalidated_snvs_file);
    }

    #-------------------------------------------------
    # look at the new calls that were found in validation, but not in the first build (Case 2 above)
    # in the case of extension experiments, with no previous build, this will be all variants
    $self->debug_message("Gathering new sites...");

    if (-s $current_snv_file){
        return 1 unless $self->gather_new_sites($current_snv_file);
    } else {
        $self->debug_message("No variants found that were called in the validation, but not found in original genomes");
    }

    ###process the validated list for examination
    $self->process_validation_list($validated_snvs_file);

    return 1;
}

sub generate_clean_file {
    my $self = shift;
    my $file = shift;

    my ($type) = $file =~ /(snvs|indels)/;
    unless($type) {
        die $self->error_message('Cannot determine variant type of file ' . $file);
    }

    my $clean_file = join('/', $self->output_dir, $type, "$type.hq.clean.bed");
    my $clean_fh = Genome::Sys->open_file_for_writing($clean_file);
    my $input_fh = Genome::Sys->open_file_for_reading($file);

    my %dups;
    while( my $line = $input_fh->getline ) {
        chomp($line);
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
        if($ref =~ /\//){ #BED file
            ( $ref, $var ) = split(/\//, $ref);
        }
        $ref =~ s/[0*]/-/g;
        $var =~ s/[0*]/-/g;

        my @vars;
        if($var !~ /-/) {
            @vars = Genome::Info::IUB->variant_alleles_for_iub($ref, $var);
        } else {
            @vars = ($var);
        }

        foreach my $v (@vars){
            unless (exists($dups{join("\t",($chr, $start, $stop, $ref, $v ))})){
                print $clean_fh join("\t",($chr, $start, $stop, $ref, $v )) . "\n";
            }
            $dups{join("\t",($chr, $start, $stop, $ref, $v ))} = 1;
        }
    }
    $clean_fh->close();
    $input_fh->close();

    return $clean_file;
}

sub parse_prior_variants {
    my $self = shift;

    my %prev_calls;
    my $variant_list;
    if(defined($self->variant_list)){
        $variant_list = $self->variant_list;
    } else { #look in model
        if(defined($self->build->snv_variant_list)){
            $variant_list = join("/", $self->build->snv_variant_list->output_dir, "snvs.hq.bed");
            $self->variant_list_is_one_based(0);
        }
    }

    if($variant_list){
        unless(-e $variant_list) {
            die $self->error_message('Could not find requested variant list: ' . $variant_list);
        }

        $self->debug_message("using variant list: $variant_list");
        my $variant_fh = Genome::Sys->open_file_for_reading($variant_list);

        while( my $line = $variant_fh->getline )
        {
            chomp($line);
            #handle either 5 col (Ref\tVar) or 4 col (Ref/Var) bed
            my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
            if($ref =~ /\//){
                ( $ref, $var ) = split(/\//, $ref);
            }
            $ref =~ s/[0*]/-/g;
            $var =~ s/[0*]/-/g;

            if(($ref eq "-") || ($var eq "-")){
                my $key = join("\t",($chr, $start, $stop, $ref, $var ));
                if($self->variant_list_is_one_based){
                    $prev_calls{annoToBed($key)} = 0;
                } else {
                    $prev_calls{$key} = 0;
                }
            } else {
                my @vars = Genome::Info::IUB->variant_alleles_for_iub($ref, $var);
                foreach my $v (@vars){
                    my $key = join("\t",($chr, $start, $stop, $ref, $v ));
                    if($self->variant_list_is_one_based){
                        $prev_calls{annoToBed($key)} = 0;
                    } else {
                        $prev_calls{$key} = 0;
                    }
                }
            }
        }
        $variant_fh->close();

    } else {
        $self->warning_message("bed file of targeted variants not found. Assuming that there were no calls for this model (and it was an extension experiment)");
    }

    return \%prev_calls;
}

sub classify_variants {
    my $self = shift;
    my $file = shift;
    my $prev_calls = shift;

    my ($type) = $file =~ /(snvs|indels)/;
    unless($type) {
        die $self->error_message('Cannot determine variant type of file ' . $file);
    }

    my $validated_fh = Genome::Sys->open_file_for_writing(join('/',$self->output_dir, $type, "$type.validated"));
    my $new_calls_fh = Genome::Sys->open_file_for_writing(join('/',$self->output_dir, $type, "$type.newcalls"));
    my $not_found_fh = Genome::Sys->open_file_for_writing(join('/',$self->output_dir, $type, "$type.notvalidated"));

    my $input_fh = Genome::Sys->open_file_for_reading($file);

    while( my $line = $input_fh->getline ) {
        chomp($line);
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
        if($ref =~ /\//){
            ( $ref, $var ) = split(/\//, $ref);
        }

        if (exists($prev_calls->{join("\t",($chr, $start, $stop, $ref, $var ))})){
            print $validated_fh bedToAnno(join("\t",($chr, $start, $stop, $ref, $var ))) . "\n";
            delete $prev_calls->{join("\t",($chr, $start, $stop, $ref, $var ))};
        } else {
            print $new_calls_fh bedToAnno(join("\t",($chr, $start, $stop, $ref, $var ))) . "\n";
        }
    }
    $input_fh->close();

    #case 3: called in original, not in validation
    foreach my $k (keys(%$prev_calls)){
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $k );

        #only print variants of the appropriate type
        next if($type eq 'snvs' xor ($ref !~ /-/ and $var !~ /-/));

        print $not_found_fh bedToAnno(join("\t",($chr, $start, $stop, $ref, $var ))) . "\n";
    }

    $validated_fh->close();
    $new_calls_fh->close();
    $not_found_fh->close();

    return 1;
}

sub prepare_roi {
    my $self = shift;
    my $feature_list = shift;

    $self->debug_message('Using regions in feature list ' . $feature_list->id);

    my $roi_bed_file = $feature_list->file_path;
    unless ( -e $roi_bed_file ){
        die $self->error_message('Feature list ' . $feature_list->id . ' has non-existent BED file <' . $roi_bed_file . '>?');
    }

    my %target_regions;
    my $roi_fh = Genome::Sys->open_file_for_reading($roi_bed_file);
    while( my $line = $roi_fh->getline ) {
        chomp($line);
        next if $line =~ /^track name/; #FIXME should select whether targets or tiles are desired--use feature-list accessor
        my ( $chr, $start, $stop, @rest) = split( /\t/, $line );
        #remove chr if present
        $chr =~ s/^chr//g;
        $target_regions{$chr}{join("\t",($start, $stop))} = 0;
    }
    $roi_fh->close;


    return \%target_regions;
}

sub region_limit_variants {
    my $self = shift;
    my $variant_file = shift;
    my $target_regions = shift;

    $self->debug_message('Filtering out off-target regions...');

    my $on_target_file = "$variant_file.on_target";
    my $on_target_fh = Genome::Sys->open_file_for_writing($on_target_file);
    my $variant_fh = Genome::Sys->open_file_for_reading($variant_file);
    #compare the variants to the targets
    VARIANT: while( my $line = $variant_fh->getline ) {
        chomp($line);
        my ( $chr, $start, $stop, @rest ) = split( /\t/, $line );

        #if we run into huge lists, this will be slow - refactor to use joinx - TODO
        my $found = 0;
        TARGET: foreach my $pos (keys(%{$target_regions->{$chr}})){
            my ($tst, $tsp) = split("\t",$pos);
            if(intersects($start,$stop,$tst,$tsp)){
                print $on_target_fh $line, "\n";
                next VARIANT;
            }
        }
    }

    $on_target_fh->close();
    $variant_fh->close();

    return $on_target_file;
}

sub tier_variant_file {
    my $self = shift;
    my $variant_file = shift;

    unless(-s $variant_file) {
        $self->debug_message('Empty variant file.  Skipping tiering.');
        return $variant_file;
    }

    my $cmd = Genome::Model::Tools::FastTier::FastTier->create(
        tier_file_location => $self->tier_file_location,
        variant_bed_file => $variant_file,
    );

    unless($cmd->execute) {
        die $self->error_message('Failed to tier variants');
    }

    return "$variant_file.tier1";
}

sub dbsnp_filter_variant_file {
    my $self = shift;
    my $variant_file = shift;

    my ($type) = $variant_file =~ /(snv|indel)s/;
    unless($type) {
        die $self->error_message('Cannot determine variant type of file ' . $variant_file);
    }
    my $variant_bed_file = Genome::Sys->create_temp_file_path;
    my $anno_to_bed_cmd = Genome::Model::Tools::Bed::Convert::AnnotationToBed->create(
        source => $variant_file,
        output => $variant_bed_file,
    );
    unless($anno_to_bed_cmd->execute()) {
        die $self->error_message('Failed to convert variant file to BED format.');
    }

    my $filter_file = $self->dbsnp_filter;

    my $novel_bed_file = Genome::sys->create_temp_file_path;
    my $cmd = Genome::Model::Tools::Joinx::Intersect->create(
        input_file_a => $variant_bed_file,
        input_file_b => $filter_file,
        output_file => Genome::Sys->create_temp_file_path, #/dev/null
        miss_a_file => $novel_bed_file,
    );
    unless($cmd->execute) {
        die $self->error_message('Failed to execute joinx');
    }

    my $pass_filter_file = "$variant_file.novel";
    my $bed_to_anno_cmd = Genome::Model::Tools::Bed::Convert::BedToAnnotation->create(

        "${type}_file" => $novel_bed_file,
        output => $pass_filter_file,
        annotator_version => 2,
    );
    unless($bed_to_anno_cmd->execute()) {
        die $self->error_message('Failed to convert joinx result to annotation format');
    }

    return $pass_filter_file;
}

sub prepare_user_filter_sites {
    my $self = shift;
    my $filter_files = shift;

    my %filter_sites;

    my @filter_files = split(",", $filter_files);
    for my $file (@filter_files) {
        unless(Genome::Sys->check_for_path_existence($file)) {
            die $self->error_message("Could not find filter file: " . $file);
        }

        my $fh = Genome::Sys->open_file_for_reading($file);
        while( my $line = $fh->getline ) {
            chomp($line);
            #handle either 5 col (Ref\tVar) or 4 col (Ref/Var) bed
            my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
            if($ref =~ /\//){
                ( $ref, $var ) = split(/\//, $ref);
            }
            $ref =~ s/0/-/g;
            $var =~ s/0/-/g;

            my @vars = Genome::Info::IUB->variant_alleles_for_iub($ref, $var);
            foreach my $v (@vars){
                $filter_sites{join("\t",($chr, $start, $stop, $ref, $v ))} = 0;
            }
        }
        $fh->close();
    }

    return \%filter_sites;
}

sub filter_variants {
    my $self = shift;
    my $filter_sites = shift;
    my $variant_file = shift;

    my $filtered_variant_file = "$variant_file.filtered";
    my $filtered_variant_fh = Genome::Sys->open_file_for_writing($filtered_variant_file);

    my $variant_fh = Genome::Sys->open_file_for_reading($variant_file);
    while( my $line = $variant_fh->getline ) {
        chomp($line);
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
        if($ref =~ /\//){
            ( $ref, $var ) = split(/\//, $ref);
        }
        unless (exists($filter_sites->{join("\t",($chr, $start, $stop, $ref, $var ))})){
            print $filtered_variant_fh $line . "\n";
        }
    }
    $variant_fh->close();
    $filtered_variant_fh->close();

    return $filtered_variant_file;
}

sub prepare_user_filter_regions {
    my $self = shift;
    my $filter_regions = shift;

    my @files = split(",", $filter_regions);
    my $filter_region_file = Genome::Sys->create_temp_file_path;

    Genome::Sys->cat(input_files => \@files, output_file => $filter_region_file);

    return $filter_region_file;
}

sub _filter_regions {
    my $self = shift;
    my $filter_regions_file = shift;
    my $variant_file = shift;

    my $filtered_variant_file = "$variant_file.filteredReg";
    my $dev_null = Genome::Sys->create_temp_file_path(); #shellcmd doesn't currently like /dev/null
    my $intersect_cmd = Genome::Model::Tools::Joinx::Intersect->create(
        input_file_a => $variant_file,
        input_file_b => $filter_regions_file,
        miss_a_file => $filtered_variant_file,
        output_file => $dev_null,
    );
    unless($intersect_cmd->execute) {
        $self->error_message("Failed to execute joinx");
        die $self->error_message;
    }

    return $filtered_variant_file;
}

sub generate_readcounts {
    my $self = shift;
    my $file = shift;
    my $filename = basename($file);

    my $build = $self->build;
    my $ref_seq_build = $self->reference_sequence_build;

    if($build->normal_sample) {
        #get readcounts from the normal bam
        my $normal_rc_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
            bam_file => $build->normal_bam,
            output_file => $self->output_dir . "/snvs/readcounts/$filename.nrm.cnt",
            variant_file => $file,
            genome_build => $ref_seq_build->full_consensus_path('fa'),
        );
        unless ($normal_rc_cmd->execute) {
            die $self->error_message("Failed to obtain normal readcounts for file $file.");
        }
    }

    my $tumor_rc_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
        bam_file => $build->tumor_bam,
        output_file => $self->output_dir . "/snvs/readcounts/$filename.tum.cnt",
        variant_file => $file,
        genome_build => $ref_seq_build->full_consensus_path('fa'),
    );
    unless ($tumor_rc_cmd->execute) {
        die $self->error_message("Failed to obtain tumor readcounts for file $file.");
    }

    return 1;
}

sub recover_low_coverage_variants {
    my $self = shift;
    my $notvalidated_snvs_file = shift;

    my $build = $self->build;

    my $som_var_model = $self->somatic_variation_model;
    unless($som_var_model) {
        $self->error_message('No somatic variation build found--cannot recover non-validated calls');
        return;
    }

    #get original bams
    my $tumor_bam_var = $som_var_model->tumor_model->last_succeeded_build->whole_rmdup_bam_file;
    unless(-e $tumor_bam_var) {
        $self->error_message("ERROR: Somatic Variation tumor bam not found");
        return;
    }
    my $normal_bam_var = $som_var_model->normal_model->last_succeeded_build->whole_rmdup_bam_file;
    unless(-e $normal_bam_var) {
        $self->error_message("ERROR: Somatic Variation normal bam not found");
        return;
    }

    my $output_dir = $self->output_dir;
    my $ref_seq_fasta = $self->reference_sequence_build->full_consensus_path('fa');

    #get readcounts from the original normal sample
    my $normal_rc2_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
        bam_file => $normal_bam_var,
        output_file => "$output_dir/snvs/readcounts/snvs.notvalidated.orig.nrm.cnt",
        variant_file => "$output_dir/snvs/snvs.notvalidated",
        genome_build => $ref_seq_fasta,
    );
    unless ($normal_rc2_cmd->execute) {
        die $self->error_message("Failed to obtain normal readcounts.");
    }

    #get readcounts from the original tumor sample
    my $tumor_rc2_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
        bam_file => $tumor_bam_var,
        output_file => "$output_dir/snvs/readcounts/snvs.notvalidated.orig.tum.cnt",
        variant_file => "$output_dir/snvs/snvs.notvalidated",
        genome_build => $ref_seq_fasta,
    );
    unless ($tumor_rc2_cmd->execute) {
        die $self->error_message("Failed to obtain tumor readcounts.");
    }

    #read in all the validation readcounts, keep only those with poor coverage
    #require 8 reads in tumor, 6 in normal, per varscan cutoffs

    my %poorly_covered_snvs;

    if($build->normal_sample) {
        my $normal_readcount_fh = Genome::Sys->open_file_for_reading( "$output_dir/snvs/readcounts/snvs.notvalidated.nrm.cnt" );
        while( my $line = $normal_readcount_fh->getline ) {
            chomp($line);
            my ( $chr, $start, $ref, $var, $refcnt, $varcnt, $vaf) = split("\t",$line);
            if(($refcnt+$varcnt) < 6){
                $poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))} = 0;
            }
        }
        close($normal_readcount_fh);
    }

    my $tumor_readcount_fh = Genome::Sys->open_file_for_reading( "$output_dir/snvs/readcounts/snvs.notvalidated.tum.cnt" );
    while( my $line = $tumor_readcount_fh->getline ) {
        chomp($line);
        my ( $chr, $start, $ref, $var, $refcnt, $varcnt, $vaf) = split("\t",$line);
        if(($refcnt+$varcnt) < 8){
            $poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))} = 0;
        }
    }
    close($tumor_readcount_fh);

    #now, go through the original readcounts, and flag any that do have good coverage for manual review
    my $original_readcounts_fh = Genome::Sys->open_file_for_reading("$output_dir/snvs/readcounts/snvs.notvalidated.orig.nrm.cnt" );
    while( my $line = $original_readcounts_fh->getline ) {
        chomp($line);
        my ( $chr, $start, $ref, $var, $refcnt, $varcnt, $vaf) = split("\t",$line);

        if(exists($poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))})){
            if(($refcnt+$varcnt) >= 20){
                $poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))} = 1;
            }
        }
    }
    close($original_readcounts_fh);

    my $poor_count = 0;
    my @poor_output;
    my $not_validated_fh = Genome::Sys->open_file_for_reading( "$output_dir/snvs/readcounts/snvs.notvalidated.orig.tum.cnt" );
    while( my $line = $not_validated_fh->getline ) {
        chomp($line);
        my ( $chr, $start, $ref, $var, $refcnt, $varcnt, $vaf) = split("\t",$line);
        if(($refcnt+$varcnt) >= 20){
            if(defined($poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))})){
                if ($poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))} == 1){
                    push(@poor_output,join("\t",( $chr, $start-1, $start, $ref, $var )));
                    $poor_count++;
                }
            }
        }
    }
    close($not_validated_fh);

    if(@poor_output){
        my $poor_coverage_fh = Genome::Sys->open_file_for_writing("$output_dir/review/poorValCoverage.bed");
        foreach my $site (@poor_output){
            print $poor_coverage_fh $site . "\n";
        }
        close($poor_coverage_fh);
    }

    my $sample_name = $build->tumor_sample->name;
    if($poor_count > 0){
        #create the xml file for this 4-way review
        my $dump_cov_xml = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
            bams => join(",",($build->normal_bam,$build->tumor_bam,$normal_bam_var,$tumor_bam_var)),
            labels => join(",",("validation normal $sample_name","validation tumor $sample_name","original normal $sample_name","original tumor $sample_name")),
            output_file => "$output_dir/review/poorValCoverage.xml",
            genome_name => $sample_name,
            review_bed_file => "$output_dir/review/poorValCoverage.bed",
            reference_name => $self->igv_reference_name,
        );
        unless ($dump_cov_xml->execute) {
            die $self->error_message("Failed to dump IGV xml for poorly covered sites.");
        }

        $self->debug_message(join("\n",
            "--------------------------------------------------------------------------------",
            "Sites for review with poor coverage in validation but good coverage in original are here:",
            "    $output_dir/review/poorValCoverage.bed",
            "IGV XML file is here:",
            "    $output_dir/review/poorValCoverage.xml"
        ));
    }
}

sub gather_new_sites {
    my $self = shift;
    my $snv_file = shift;
    my $build = $self->build;
    my $ref_seq_fasta = $self->reference_sequence_build->full_consensus_path('fa');

    #can't run UHC if tumor-only:
    if($build->normal_sample){
        $self->debug_message("Running UHC filter...");
        #run the uhc filter to remove solid calls
        my $uhc_cmd = Genome::Model::Tools::Somatic::UltraHighConfidence->create(
            normal_bam_file => $build->normal_bam,
            tumor_bam_file => $build->tumor_bam,
            output_file => "$snv_file.passuhc",
            variant_file => "$snv_file",
            reference => $ref_seq_fasta,
            filtered_file => "$snv_file.failuhc",
        );
        unless ($uhc_cmd->execute) {
            die $self->error_message("Failed to run UHC filter.");
        }
    }


    #now get the files together for review
    $self->debug_message("Generating Review files...");
    my $revfile;
    if ( -s "$snv_file.passuhc"){
        $revfile = "$snv_file.passuhc";
    } else {
        $revfile = "$snv_file";
    }

    my $output_dir = $self->output_dir;
    my $newcalls_bed_file = "$revfile.bed";
    my $anno_to_bed_cmd = Genome::Model::Tools::Bed::Convert::AnnotationToBed->create(
        source => $revfile,
        output => $newcalls_bed_file,
    );
    unless($anno_to_bed_cmd->execute) {
        die $self->error_message('Failed to convert review file to BED format.');
    }

    my $tier1_bed_file;
    eval{$tier1_bed_file = $self->tier_variant_file($newcalls_bed_file)};
    return if $@;

    my $tier1_review_bed_file = "$output_dir/review/newcalls.bed";
    Genome::Sys->copy_file($tier1_bed_file, $tier1_review_bed_file);

    my $sample_name = $build->tumor_sample->name;

    my $bam_files;
    my $labels;
    unless($build->normal_sample){
        if(defined($self->somatic_variation_model)){
            #add tumor and normal from somatic-variation model
            my $som_var_model = $self->somatic_variation_model;
            my $tbam = $som_var_model->tumor_model->last_succeeded_build->whole_rmdup_bam_file;
            my $nbam = $som_var_model->normal_model->last_succeeded_build->whole_rmdup_bam_file;
            $bam_files = join(",",($build->tumor_bam,$tbam,$nbam));
            $labels = join(",",("validation tumor $sample_name","original tumor $sample_name","original normal $sample_name"));
        } else {
            $bam_files = join(",",($build->tumor_bam));
            $labels = join(",",("validation tumor $sample_name"));
        }
    } else {
        $bam_files = join(",",($build->normal_bam,$build->tumor_bam));
        $labels = join(",",("validation normal $sample_name","validation tumor $sample_name"));
    }


    #create the xml file for review
    my $dump_xml = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
        bams => $bam_files,
        labels => $labels,
        output_file => "$output_dir/review/newcalls.xml",
        genome_name => $sample_name,
        review_bed_file => $tier1_review_bed_file,
        reference_name => $self->igv_reference_name,
    );
    unless ($dump_xml->execute) {
        die "Failed to dump IGV xml for poorly covered sites.\n";
    }

    $self->debug_message(join("\n",
        "--------------------------------------------------------------------------------",
        "Sites to review that were not found original genomes, but were found in validation are here:",
        "    $output_dir/review/newcalls.bed.tier1",
        "IGV XML file is here:",
        "    $output_dir/review/newcalls.xml",
    ));

    return 1;
}

sub process_validation_list {
    my $self = shift;
    my $variant_file = shift;

    my $variant_bed_file = "$variant_file.bed";
    my $anno_to_bed_cmd = Genome::Model::Tools::Bed::Convert::AnnotationToBed->create(
        source => $variant_file,
        output => $variant_bed_file,
    );
    unless($anno_to_bed_cmd->execute()) {
        die $self->error_message('Failed to convert variant file back to BED format.');
    }

    my $tier1_bed = $self->tier_variant_file($variant_bed_file);

    if(-s $tier1_bed) { #only annotate if we have something to annotate!
        my $annotation_build = $self->build->annotation_build;
        my $annotator_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
            use_version => 2,
            variant_bed_file => $tier1_bed,
            annotation_filter => 'top',
            build_id => $annotation_build->id,
            output_file => "$tier1_bed.anno"
        );
        unless($annotator_cmd->execute()) {
            die $self->error_message('Failed to annotated variants.');
        }
    }

    return 1;
}

1;
