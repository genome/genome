package Genome::Model::MutationalSignificance::Command::CreateMafFile;

use strict;
use warnings;

use Genome;
use Sort::Naturally qw(nsort);

class Genome::Model::MutationalSignificance::Command::CreateMafFile {
    is => ['Command::V2'],
    has_input => [
        somatic_variation_build => {
            is => 'Genome::Model::Build::SomaticVariation',
        },
        output_dir => {
            is => 'Text',
        },
        include_skipped_cosmic_dbsnps => {
            is => 'Boolean',
            default_value => 1,
        },
        include_regulatory => {
            is => 'Boolean',
            default_value => 0,
        },
        include_ensembl_annot => {
            is => 'Boolean',
            default_value => 1,
        },
        filter_on_regulomedb => {
            is => 'Boolean',
            default_value => 1,
        },
        cosmic_dir => {
            is => 'Path',
            doc => 'cosmic amino acid mutation database folder',
            default => Genome::Sys->dbpath('cosmic','latest'),
        },
        exclude_pindel_only_indels => {
            is => 'Boolean',
            default_value => 0,
        },
        review_file_dir => {
            is => 'UR::Value::DirectoryPath',
            doc => 'Path to directory of variant files with reviews.  Any variant with a review status other than S or V will be ignored.',
            is_optional => 1,
        },
        use_tier_1 => {
            is => 'Boolean',
            doc => 'Include tier 1 in the analysis',
            default_value => 1,
        },
        use_tier_2 => {
            is => 'Boolean',
            doc => 'Include tier 2 in the analysis',
            default_value => 0,
        },
        use_tier_3 => {
            is => 'Boolean',
            doc => 'Include tier 3 in the analysis',
            default_value => 0,
        },
        use_tier_4 => {
            is => 'Boolean',
            doc => 'Include tier 4 in the analysis',
            default_value => 0,
        },
        regulatory_columns_to_check => {
            is => 'String',
            doc => 'Names of columns in the regulatory annotation file to check',
            is_optional => 1,
            is_many => 1,
        },
        use_indels => {
            is => 'Boolean',
            doc => 'Whether to use indels',
            default => 1,
        },
    ],
    has_output => [
        maf_file => {},
    ],
};

sub execute {
    my $self = shift;

    $self->status_message("CreateMafFile for build ".$self->somatic_variation_build->id);
    #my $snv_file = $self->output_dir."/".$self->somatic_variation_build->id.".snv.anno";

    #Exclude reviewed variants that didn't pass review
    my %failed_variants;
    my $id = $self->somatic_variation_build->tumor_build->model->subject->extraction_label;
    my @review_files;
    if ($self->review_file_dir) {
        @review_files = glob($self->review_file_dir."/".$id."*.reviewed.csv");
    }
    foreach my $file (@review_files) {
        my @review_lines = `cat $file` if (-e $file);
        chomp @review_lines;
        foreach my $line (@review_lines){
            my ($chr, $start, $stop, undef, undef, $call) = split(/\t/, $line);
            next if($chr eq 'Chr'); #skip header
            unless (!defined $call or $call =~ m/^s*$/ or $call =~ m/[SVsv]/){
                $failed_variants{"$chr\t$start\t$stop"} = $line;
            }
        }
    }

    #Deduplicate and sort the snv file (copied from gmt capture manual-review)
    my ($snv_anno_fh, $snv_anno_file) = Genome::Sys->create_temp_file;

    my @tiers_to_use;
    if ($self->use_tier_1) {
        push @tiers_to_use, 1;
    }
    if ($self->use_tier_2) {
        push @tiers_to_use, 2;
    }
    if ($self->use_tier_3) {
        push @tiers_to_use, 3;
    }
    if ($self->use_tier_4) {
        push @tiers_to_use, 4;
    }
    my $version = 1;
    my $snv_anno_with_rsid = $self->output_dir."/rsid";
    foreach my $tier (@tiers_to_use){
        my $snv_anno_top = $self->somatic_variation_build->data_set_path("effects/snvs.hq.tier$tier",$version,"annotated.top.header");
        my $snv_anno_rsid = $self->somatic_variation_build->data_set_path("effects/snvs.hq.tier$tier.rsid", $version, "annotated.top");
        `cat $snv_anno_rsid >> $snv_anno_with_rsid`;
        my $snv_regulatory = $self->somatic_variation_build->data_set_path("effects/snvs.hq.tier$tier",$version, "annotated.top.header.regulatory");
        my $snv_anno = $snv_anno_top;
        if ($self->include_regulatory and -s $snv_regulatory) {

            $snv_anno = $self->output_dir."/".$self->somatic_variation_build->id."tier$tier.merged_anno";
            my %params = (
                tgi_anno_file => $snv_anno_top,
                regulatory_file => $snv_regulatory,
                regulatory_columns_to_check => [$self->regulatory_columns_to_check],
                output_file => $snv_anno,
                annotation_build => $self->somatic_variation_build->annotation_build,
                include_ensembl_annot => $self->include_ensembl_annot,
            );

            if ($self->filter_on_regulomedb) {
                my $rdb_file = $self->somatic_variation_build->data_set_path("effects/snvs.hq.regulomedb", 1, "full");
                unless (-s $rdb_file) {
                    $self->error_message("No regulomedb file detected");
                    return;
                }

                $params{regulome_db_file} = $rdb_file;
            }

            my $rv = Genome::Model::MutationalSignificance::Command::MergeAnnotations->execute(
                %params
            );
           
            unless ($rv) {
                $self->error_message("Failed to merge annotations for tier $tier");
                return $rv;
            }
        }

        my @snv_lines = `cat $snv_anno`;
        chomp @snv_lines;

        #Remove any variants that failed review
        for my $line ( @snv_lines )
        {
            my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
            if (%failed_variants and defined $failed_variants{"$chr\t$start\t$stop"}) {
                next;
            }
            else {
                $snv_anno_fh->print($line."\n");
            }
        }
    }
    $snv_anno_fh->close;

    $snv_anno_file = $self->sort_and_write_to_temp_file($snv_anno_file);
    my $indel_anno_file;
    if ($self->use_indels) {
        #Exclude pindel-only indels
        my %uniq_to_pindel;
        if ($self->exclude_pindel_only_indels) {
            my ($gatk_indels) = glob($self->somatic_variation_build->data_directory."/variants/indel/gatk-somatic-indel-*/indels.hq.bed");
            my ($varscan_indels) = glob($self->somatic_variation_build->data_directory."/variants/indel/varscan-somatic-*/varscan-high-confidence-indel-*/indels.hq.bed");
            my ($pindel_indels) = glob($self->somatic_variation_build->data_directory."/variants/indel/pindel-*/pindel-somatic-calls-*/pindel-vaf-filter-*/pindel-read-support-*/indels.hq.bed");

            my %gatk_varscan_lines = map {chomp; s/\*/0/; $_=>1} `cut -f 1-4 $varscan_indels $gatk_indels`;
            %uniq_to_pindel = map {chomp; $_=>1} `cut -f 1-4 $pindel_indels` if( -e $pindel_indels );
            foreach my $indel ( keys %gatk_varscan_lines ) {
                delete $uniq_to_pindel{$indel} if( defined $uniq_to_pindel{$indel} );
            }
        }

        my @indel_lines;
        foreach my $tier (@tiers_to_use){
            my $indel_anno = $self->somatic_variation_build->data_set_path("effects/indels.hq.tier$tier", $version, "annotated.top");

            my @tier_indel_lines = `cat $indel_anno`;
            chomp @tier_indel_lines;
            @indel_lines = (@indel_lines, @tier_indel_lines);
        }

        #Deduplicate indel files and remove pindel-only and failed variants if needed
        my %review_lines;
        my $indel_cnt = 0;
        for my $line (@indel_lines) {
            my ($chr, $start, $stop, $ref, $var) = split(/\t/, $line);
            my ($base0_start, $base0_stop) = ($start-1,$stop-1);
            my $refvar = ($ref eq '-' ? "0/$var": "$ref/0");

            if (%failed_variants and 
                (defined $failed_variants{"$chr\t$start\t$stop"} or
                    defined $failed_variants{"$chr\t$base0_start\t$base0_stop"})) {
                next;
            }
            elsif ($self->exclude_pindel_only_indels and (defined $uniq_to_pindel{"$chr\t$base0_start\t$stop\t$refvar"} or 
                    defined $uniq_to_pindel{"$chr\t$start\t$base0_stop\t$refvar"} )){
                $review_lines{pindels}{$chr}{$start}{$stop} = $line;
            }
            else {
                ++$indel_cnt unless( defined $review_lines{indels}{$chr}{$start}{$stop} );
                $review_lines{indels}{$chr}{$start}{$stop} = $line;
            }
        }

        #Write indel lines
        my ($indel_anno_fh, $indel_anno_file_temp) = Genome::Sys->create_temp_file;
        foreach my $chr (nsort keys %{$review_lines{indels}}) {
            foreach my $start (sort {$a <=> $b} keys %{$review_lines{indels}{$chr}}) {
                foreach my $stop (sort {$a <=> $b} keys %{$review_lines{indels}{$chr}{$start}}) {
                    $indel_anno_fh->print($review_lines{indels}{$chr}{$start}{$stop}."\n");
                }
            }
        }
        $indel_anno_fh->close;
        $indel_anno_file = $indel_anno_file_temp;
    }

    #TODO: Option to remove any that fail manual review
    #TODO: Check count of variants to review and set aside if too many (maybe put this in merge-maf?)
    
    #Include the ultra-high-confidence snvs
    #my $uhc_cmd = Genome::Model::Tools::Somatic::UltraHighConfidence->create(
    #    normal_bam_file => $self->somatic_variation_build->normal_bam,
    #    tumor_bam_file => $self->somatic_variation_build->tumor_bam,
    #    variant_file => $snv_anno_file,
    #    output_file => $snv_file,
    #    filtered_file => $self->output_dir."/".$self->somatic_variation_build->id.".not_uhc.anno",
    #    reference => $self->somatic_variation_build->reference_sequence_build->fasta_file,
    #);

    #my $uhc_result = $uhc_cmd->execute;

    #Add in the dbSnp variants that appear in COSMIC
    if ($self->include_skipped_cosmic_dbsnps) {
        my $known_cosmic_sites = Genome::Sys->create_temp_file_path;
        my $cosmic_db = $self->cosmic_dir."/Cosmic_Database.tsv";
        my $parse_cosmic_cmd = "cut -f 7,14-16 $cosmic_db | egrep \"Confirmed|Previously Observed\" | cut -f 2-4 | perl -ane '\$F[0]=~s/^23\$/X/; \$F[0]=~s/^24\$/Y/; if(\$F[0]=~m/^\\w+\$/ && \$F[1]=~m/^\\d+\$/ && \$F[2]=~m/^\\d+\$/) {print join(\"\t\",\@F).\"\n\"}' | joinx sort -s > $known_cosmic_sites";
        Genome::Sys->shellcmd(cmd => $parse_cosmic_cmd);

        my %cosmic_loci = map {chomp; ($_,1) } `cat $known_cosmic_sites`;

        my @matching_skipped_loci;
        my $skipped_bed = $self->somatic_variation_build->data_set_path("effects/snvs.hq.previously_detected.tier1",2,"bed");
        next unless (-e $skipped_bed); #nothing was skipped
        my $skipped = Genome::Sys->create_temp_file_path;
        my $annotation_build = $self->somatic_variation_build->annotation_build;
        unless ($annotation_build) {
            print "annotation build was null for somatic variation build ".$self->somatic_variation_build->id."\n";
        }
        my $annotation_build_id = $self->somatic_variation_build->annotation_build->id; 
        my $rv = `gmt annotate transcript-variants --variant-bed-file $skipped_bed --output-file $skipped --build-id $annotation_build_id --annotation-filter top --accept-reference-IUB-codes`;
        my @skipped_loci = `cat $skipped`;

        foreach my $locus (@skipped_loci) {
            my ($chr, $start, $stop) = split(/\t/, $locus);
            if (defined $cosmic_loci{"$chr\t$start\t$stop"}) {
                push (@matching_skipped_loci, $locus);
            }
        }
        my $fh = IO::File->new( $snv_anno_file, "a");
        $fh->print(@matching_skipped_loci);
        $fh->close;

        $snv_anno_file = $self->sort_and_write_to_temp_file($snv_anno_file);
    }
    
    my %params = (
        snv_file => $snv_anno_file,
        snv_annotation_file => $snv_anno_file,
        snv_anno_file_with_rsid => $snv_anno_with_rsid,
        genome_build => $self->somatic_variation_build->reference_sequence_build->version,
        tumor_sample => $self->somatic_variation_build->tumor_build->model->subject->extraction_label, #TODO verify
        normal_sample => $self->somatic_variation_build->normal_build->model->subject->extraction_label, #TODO verify
        output_file => $self->output_dir."/".$self->somatic_variation_build->id.".maf",
    );
    if ($self->use_indels) {
        $params{indel_file} = $indel_anno_file;
        $params{indel_annotation_file} = $indel_anno_file;
    }
    my $create_maf_cmd = Genome::Model::Tools::Capture::CreateMafFile->create(%params);

    my $create_maf_result = $create_maf_cmd->execute;

    $self->maf_file($self->output_dir."/".$self->somatic_variation_build->id.".maf");
    my $status = "Created MAF file ".$self->maf_file;
    $self->status_message($status);
    return 1;
}

sub sort_and_write_to_temp_file {
    my ($self,$file) = @_;

    my @snv_lines = `cat $file`;

    chomp @snv_lines;
    # Store the variants into a hash to help sort variants by loci, and to remove duplicates
    my %review_lines = ();
    for my $line ( @snv_lines )
    {
        my ( $chr, $start, $stop, $ref, $var, $type, $call ) = split( /\t/, $line );
        $review_lines{snvs}{$chr}{$start}{$stop}{$call} = $line;
    }
    my ($snv_anno_fh, $snv_anno_file) = Genome::Sys->create_temp_file;
    for my $chr ( nsort keys %{$review_lines{snvs}} )
    {
        for my $start ( sort {$a <=> $b} keys %{$review_lines{snvs}{$chr}} )
        {
            for my $stop ( sort {$a <=> $b} keys %{$review_lines{snvs}{$chr}{$start}} )
            {
                for my $call (keys %{$review_lines{snvs}{$chr}{$start}{$stop}}) {
                    $snv_anno_fh->print( $review_lines{snvs}{$chr}{$start}{$stop}{$call}, "\n" );
                }
            }
        }
    }
    $snv_anno_fh->close;
    return $snv_anno_file;
}

1;
