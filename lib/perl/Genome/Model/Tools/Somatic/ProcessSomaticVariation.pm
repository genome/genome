package Genome::Model::Tools::Somatic::ProcessSomaticVariation;

use warnings;
use strict;
use IO::File;
use Genome;
use Sort::Naturally qw(nsort);
use Genome::Info::IUB;
use Spreadsheet::WriteExcel;
use File::Slurp qw(read_dir);
use File::Basename qw(basename);
use Carp qw(confess);
use Params::Validate;
use YAML;

class Genome::Model::Tools::Somatic::ProcessSomaticVariation {
    is => 'Command::V2',
    has_input => [
        somatic_variation_build => {
            is => 'Genome::Model::Build::SomaticVariation',
            doc => 'Somactic Variation build',
        },
        igv_reference_name =>{
            is => 'Text',
            doc => "name of the igv reference to use for review",
            example_values => ["reference_build36","b37","mm9"],
        },
    ],
    has_optional_input => [
        # make pp option
        restrict_to_target_regions =>{
            is => 'Boolean',
            default => 1,
            doc => "only keep snv calls within the target regions. These are pulled from the build if possible",
        },
        # make pp option
        target_regions =>{
            is => 'String',
            doc => "path to a target file region. Used in conjunction with --restrict-to-target-regions to limit sites to those appearing in these regions",
        },
        # pp option
        required_snv_callers => {
            is => 'Number',
            doc => "Number of independent algorithms that must call a SNV. If set to 1 (default), all calls are used",
            default => 1,
        },
        # pp option
        tiers_to_review => {
            is => 'String',
            doc => "comma-separated list of tiers to include in review. (e.g. 1,2,3 will create bed files with tier1, tier2, and tier3)",
            default => 1,
        },
        # goes away
        sample_name =>{
            is => 'Text',
            is_mutable => 0,
            calculate => q{ return $self->somatic_variation_build->get_subject_name; },
            doc => "override the sample name on the build and use this name instead",
        }
    ],
    has_optional_output => [
        report => {
            is => 'Path',
            calculate =>  q{ File::Spec->join($_output_dir, 'snvs.indels.annotated') },
            calculate_from => ['_output_dir'],
        },
        report_xls => {
            is => 'Path',
            calculate =>  q{ $report . '.xls' },
            calculate_from => ['report'],
        },
        review_dir => {
            is => 'Path',
            calculate => q{ File::Spec->join($_output_dir, 'review') },
            calculate_from => ['_output_dir'],
        },
        review_bed => {
            is => 'Path',
            calculate =>  q{ File::Spec->join($review_dir, $sample_name . '.bed') },
            calculate_from => ['review_dir', 'sample_name' ],
        },
        review_xml => {
            is => 'Path',
            calculate =>  q{ File::Spec->join($review_dir, $sample_name . '.xml') },
            calculate_from => ['review_dir', 'sample_name' ],
        },
    ],
    has_transient_optional => [
        _build_dir => {
            is => 'Text',
            is_mutable => 0,
            calculate => q{ $self->somatic_variation_build->data_directory },
        },
        _output_dir => {
            is => 'Text',
            doc => "Directory where output will be stored",
        },
        _allocation => {
            is => 'Genome::Disk::Allocation',
        },
    ],
};


sub help_detail {
  return <<HELP;
Given a SomaticVariation build, this tool will gather the resulting variants, remove
off-target sites, tier the variants, optionally filter them, etc. Calls are prepped for
manual review in the review/ directory.
HELP
}

sub _doc_authors {
  return <<AUTHS;
 Chris Miller
 Susanna Siebert
 David Morton
AUTHS
}


sub execute {
    my $self = shift;

    $self->status_message("Processing model with sample_name: " . $self->sample_name);

    $self->_output_dir($self->create_allocation());

    # Check if the necessary files exist in this build and put them in the processing location
    my $snv_file   = $self->stage_snv_file();
    my $indel_file = $self->stage_indel_file();

    #--------------------------------------------------------------
    # munge through SNV file to remove duplicates and fix IUB codes
    $snv_file = $self->cleanFile($snv_file, $self->snvs_dir);
    $indel_file = $self->cleanFile($indel_file, $self->indels_dir);

    #-------------------------------------------------------
    # remove regions called by less than the required number of callers
    $snv_file = $self->removeUnsupportedSites($snv_file);

    #-------------------------------------------------
    # filter out the off-target regions, if target regions are available
    if ($self->restrict_to_target_regions) {
        $snv_file = $self->_filter_off_target_regions($snv_file, $self->snvs_dir);
        $indel_file = $self->_filter_off_target_regions($indel_file, $self->indels_dir);
    }

    #------------------------------------------------------
    # do annotation
    $snv_file   = $self->doAnnotation($snv_file, $self->snvs_dir);
    $indel_file = $self->doAnnotation($indel_file, $self->indels_dir);

    #-------------------------------------------------------
    # add tiers
    $self->status_message("Adding tiers");
    $snv_file   = $self->addTiering($snv_file);
    $indel_file = $self->addTiering($indel_file);

    #convert back to annotation format (1-based)
    $snv_file = bedFileToAnnoFile($snv_file);
    $indel_file = bedFileToAnnoFile($indel_file);

    #----------------------------------------------------
    # add dbsnp/gmaf
    ($snv_file, $indel_file) = $self->_add_dbsnp_and_gmaf($snv_file, $indel_file);

    $self->status_message("Getting read counts");
    $snv_file   = $self->read_counts($snv_file, $self->snvs_dir);
    $indel_file = $self->read_counts($indel_file, $self->indels_dir);

    #------------------------------------------------------
    # combine the files into one master table
    $self->_create_master_files($snv_file, $indel_file);

    #------------------------------------------------------
    # now get the files together for review
    $self->_create_review_files();

    $self->reallocate();
    #$self->create_symlinks_in_build_dir();

    return 1;
}

sub create_allocation {
    my $self = shift;

    my $build = $self->somatic_variation_build;
    my %allocation_params = (
        disk_group_name => $ENV{GENOME_DISK_GROUP_MODELS},
        allocation_path => File::Spec->join('model_data', 'ProcessSomaticVariation', $build->id),
        kilobytes_requested => 4_000_000,
        owner_class_name => $build->class,
        owner_id => $build->id,
    );

    my $allocation = Genome::Disk::Allocation->allocate(
        %allocation_params);
    unless ($allocation) {
        confess $self->error_message(
                "Failed to get disk allocation with params:\n" .
                YAML::Dump(%allocation_params));
    }

    my $absolute_path = $allocation->absolute_path;
    unless (-d $absolute_path) {
        $allocation->delete;
        confess $self->error_message(
                "Path $absolute_path doesn't exist!");
    }
    $self->_allocation($allocation);

    return $absolute_path;
}

sub reallocate {
    my $self = shift;

    $self->_allocation->reallocate();
}

sub bedFileToAnnoFile{
    my ($file, $outfile) = @_;

    #remove bed from name
    my $newfile = $file;
    $newfile =~ s/\.bed//g;

    if ($outfile) {
        $newfile = $outfile;
    }

    my $outFh = Genome::Sys->open_file_for_writing($newfile);
    my $inFh  = Genome::Sys->open_file_for_reading($file);
    while (my $line = $inFh->getline) {
        chomp($line);
        if ($line =~ /^chrom/) {
            print $outFh $line . "\n";
            next;
        }
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        if ($ref =~ /\//) {
            my @alleles = split("/", $ref);
            print $outFh bedToAnno(join("\t",($chr, $start, $stop, $alleles[0], $alleles[1]))) . "\t" . join("\t",($var,@rest)) . "\n";
        }
        else {
            print $outFh bedToAnno(join("\t",($chr, $start, $stop, $ref, $var))) . "\t" . join("\t",@rest) . "\n";
        }
    }
    close($outFh);
    close($inFh);
    return $newfile;
}

sub annoFileToBedFile{
    my ($file, $outfile) = @_;

    #add bed to name
    my $newfile = $file . ".bed";
    if ($outfile) {
        $newfile = $outfile;
    }

    my $outFh = Genome::Sys->open_file_for_writing($newfile);
    my $inFh  = Genome::Sys->open_file_for_reading($file);
    while (my $line = $inFh->getline) {
        chomp($line);
        if ($line =~ /^chrom/) {
            print $outFh $line . "\n";
            next;
        }
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        print $outFh annoToBed(join("\t",($chr, $start, $stop, $ref, $var))) . "\t" . join("\t",@rest) . "\n";
    }
    close($outFh);
    close($inFh);
    return $newfile;
}

sub annoFileToSlashedBedFile{
    my ($file, $outfile) = @_;

    my $newfile = $file . ".bed";
    if ($outfile) {
        $newfile = $outfile;
    }

    my $outFh = Genome::Sys->open_file_for_writing($newfile);
    my $inFh  = Genome::Sys->open_file_for_reading($file);
    while (my $line = $inFh->getline) {
        chomp($line);
        if ($line =~ /^chrom/) {
            next;
        }
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        my $bed = annoToBed(join("\t",($chr, $start, $stop, $ref, $var)));
        my @bedline = split(/\t/, $bed);
        print $outFh join("\t",(@bedline[0..2],"$bedline[3]/$bedline[4]")) . "\n";
    }
    close($outFh);
    close($inFh);
    return $newfile;
}

sub bedToAnno{
    my ($chr, $start, $stop, $ref, $var) = split("\t", $_[0]);
    #$self->status_message(join("|",($chr, $start, $stop, $ref, $var)));
    if ($ref =~ /^[-0*]/) { #indel INS
        $stop++;
    }
    else { #indel DEL or SNV
        $start++;
    }
    return join("\t", ($chr, $start, $stop, $ref, $var));
}


sub annoToBed{
    my ($chr, $start, $stop, $ref, $var) = split("\t", $_[0]);
    if ($ref =~ /^[-*0]/) { #indel INS
        $stop--;
    }
    else { #indel DEL or SNV
        $start--;
    }
    #handle 5 col or 4 col ref/var
    if (defined($var)) {
        return join("\t", ($chr, $start, $stop, $ref, $var));
    }
    else {
        return join("\t", ($chr, $start, $stop, $ref));
    }
}

sub intersects{
    my ($st, $sp, $st2, $sp2) = @_;
    if ((($sp2 >= $st) && ($sp2 <= $sp)) ||
        (($sp >= $st2) && ($sp <= $sp2))) {
        return 1;
    }
    return 0;
}

sub fixIUB{
    my ($ref, $var) = @_;
    my @vars = Genome::Info::IUB->variant_alleles_for_iub($ref, $var);
    return @vars;
}

sub addName{
    my ($file, $name) = @_;
    if ($file =~ /\.bed$/) {
        $file =~ s/\.bed//g;
        return $file . "." . $name . ".bed";
    }
    else {
        return $file . "." . $name;
    }
}

#read in the file, output a cleaned-up version
sub cleanFile {
    my ($self, $file, $directory) = @_;

    my $newfile = $self->result_file_path(
        input_file_path => $file,
        suffix => 'clean',
        directory => $directory,
    );

    my %dups;

    my $outFh = Genome::Sys->open_file_for_writing($newfile);
    my $inFh  = Genome::Sys->open_file_for_reading($file);
    while (my $line = $inFh->getline) {
        chomp($line);
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        if ($ref =~ /\//) {
            ( $ref, $var ) = split(/\//, $ref);
        }

        $ref =~ s/0/-/g;
        $var =~ s/0/-/g;
        $ref =~ s/\*/-/g;
        $var =~ s/\*/-/g;

        my @vars = ($var);
        unless ($ref =~ /-/ || $var =~ /-/) { #fixiub doesn't handle indels
            @vars = fixIUB($ref, $var);
        }

        foreach my $v (@vars) {
            unless (exists($dups{join("\t",($chr, $start, $stop, $ref, $v ))})) {
                print $outFh join("\t",($chr, $start, $stop, $ref, $v )) . "\n";
            }
            $dups{join("\t",($chr, $start, $stop, $ref, $v ))} = 1;
        }
    }
    close($outFh);
    close($inFh);
    Genome::Sys->shellcmd( cmd => "joinx sort -i $newfile >$newfile.tmp" );
    Genome::Sys->shellcmd( cmd => "mv -f $newfile.tmp $newfile");
    return $newfile;
}


sub getSiteHash  {
    my $self = shift;
    my $filter_string = shift;

    my @filters = split(",", $filter_string);
    my %filterSites;

    foreach my $filterfile (@filters) {
        if (-s $filterfile) {
            #store sites to filter out in a hash
            my $inFh = Genome::Sys->open_file_for_reading($filterfile);
            while (my $line = $inFh->getline) {
                chomp($line);
                #handle either 5 col (Ref\tVar) or 4 col (Ref/Var) bed
                my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
                if ($ref =~ /\//) {
                    ( $ref, $var ) = split(/\//, $ref);
                }
                $ref =~ s/0/-/g;
                $var =~ s/0/-/g;
                $ref =~ s/\*/-/g;
                $var =~ s/\*/-/g;

                my @vars = fixIUB($ref, $var);
                foreach my $v (@vars) {
                    $filterSites{join("\t",($chr, $start, $stop, $ref, $v ))} = 0;
                }
            }
            close($inFh);
        }
        else {
            confess $self->error_message("filter sites file doesn't exist or has zero size: " . $filterfile);
        }
    }
    return \%filterSites;
}

##TODO - these two functions are brittle. There needs to be a better way to grab calls for specific callers. Ideally, from the VCF...
sub getVcfFile{
    my $self      = shift;
    my $dir       = shift;

    my $prefix = $self->_build_dir . "/variants/snv/";

    if ($dir=~/strelka/) {
        if (-s "$prefix/$dir/snvs.vcf.gz") {
            return "$prefix/$dir/snvs.vcf.gz";
        }

    }
    elsif ($dir=~/varscan/) {
        if (-s "$prefix/$dir/varscan-high-confidence-v1-d41d8cd98f00b204e9800998ecf8427e/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/snvs.vcf.gz") {
            return "$prefix/$dir/varscan-high-confidence-v1-d41d8cd98f00b204e9800998ecf8427e/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/snvs.vcf.gz";
        }
    }
    elsif ($dir=~/sniper/) {
        if (-s "$prefix/$dir/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/somatic-score-mapping-quality-v1-39b60f48b6f8c9e63436a5424305e9fd/snvs.vcf.gz") {
            return "$prefix/$dir/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/somatic-score-mapping-quality-v1-39b60f48b6f8c9e63436a5424305e9fd/snvs.vcf.gz";
        }
    }
    else {
        confess $self->error_message("Don't know how to find the calls for $prefix/$dir. $dir is neither 'strelka', 'varscan' or 'sniper'.");
    }
    confess $self->error_message("Couldn't find the snvs.vcf.gz file under directory $dir");
}


sub removeUnsupportedSites {
    my ($self, $snv_file) = @_;

    if ($self->required_snv_callers == 1) {
        return $snv_file;
    }
    $self->status_message("Removing sites supported by less than %s callers", $self->required_snv_callers);

    #hash all of the sites
    my $sites = $self->getSiteHash($snv_file);

    #Look for the callers
    my @dirs = map {basename($_) } glob(
        File::Spec->join($self->_build_dir, 'variants', 'snv', '*'),
    );
    #remove non-caller dirs
    @dirs = grep{ ! /^intersect|^union|^samtools/ } @dirs;

    #count the number of callers that called each site from the vcfs
    for my $dir (@dirs) {
        my $file = $self->getVcfFile($dir);
        my $ifh = Genome::Sys->open_gzip_file_for_reading($file);

        while (my $line = $ifh->getline) {
            chomp $line;
            next if $line =~ /^#/;
            my ($chr, $pos, $id, $ref, $var, @rest) = split /\t/, $line;
            my @vars = split(",", $var);
            for my $v (@vars) {
                my $key = join("\t",($chr, $pos-1, $pos, $ref, $v));
                if (defined($sites->{$key})) {
                    $sites->{$key}++;
                }
            }
        }
        $ifh->close;
    }

    my $result_file_path = $self->result_file_path(
        input_file_path => $snv_file,
        suffix          => "gt" . $self->required_snv_callers . "callers",
        directory    => $self->snvs_dir,
    );

    my $ofh = Genome::Sys->open_file_for_writing($result_file_path);
    #read the snv_file again to preserve order, traiing fields, etc.
    my $ifh = Genome::Sys->open_file_for_reading($snv_file);
    while (my $line = $ifh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $ref, $var, @rest) = split /\t/, $line;
        my $key = join("\t",($chr, $start, $stop, $ref, $var));

        if (!defined($sites->{$key})) {
            $self->status_message("wut?: " . $key);
        }
        if ($sites->{$key} >= $self->required_snv_callers) {
            print $ofh $line . "\n";
        }
    }
    close($ifh);
    close($ofh);
    return $result_file_path;
}

sub result_file_path {
    my $self   = shift;
    my %params = Params::Validate::validate(@_,
        {
            input_file_path => Params::Validate::SCALAR,
            suffix => Params::Validate::SCALAR,
            directory => Params::Validate::SCALAR,
        },
    );

    my $input_file_name  = basename($params{input_file_path});
    my $output_file_name = addName($input_file_name, $params{suffix});
    my $result_file_path = File::Spec->join($params{directory}, $output_file_name);

    return $result_file_path;
}

sub doAnnotation {
    my $self         = shift;
    my $file         = shift;
    my $directory    = shift;

    my $annotation_build_name = $self->annotation_build->name;

    if ($file =~ /.bed$/) {
        $file = bedFileToAnnoFile($file);
    }

    my $newfile = $self->result_file_path(
        input_file_path => $file,
        suffix          => "anno",
        directory       => $directory,
    );

    #handle zero size files
    if (-z $file ) {
        Genome::Sys->shellcmd( cmd => "touch $newfile" );
        return $newfile;
    }

    my $anno_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
        variant_file          => $file,
        output_file           => $newfile,
        reference_transcripts => $annotation_build_name,
        annotation_filter     => "top",
    );
    unless ($anno_cmd->execute) {
        confess $self->error_message("Failed to annotate variants successfully.");
    }
    return $newfile;
}


sub addTiering {
    my $self               = shift;
    my $file               = shift;
    my $outfile            = shift;

    my $tier_file_location = File::Spec->join($self->annotation_build->data_directory, 'annotation_data', 'tiering_bed_files_v3');

    unless ($file =~ /\.bed/) {
        $file = annoFileToBedFile($file);
    }

    my $newfile;
    if ($outfile) {
        $newfile = $outfile;
    }
    else {
        $newfile = addName($file, "tiered");
    }

    #handle zero size files
    if (-z $file) {
        Genome::Sys->shellcmd( cmd => "touch $newfile" );
        return $newfile;
    }

    my $tier_cmd = Genome::Model::Tools::FastTier::AddTiers->create(
        input_file         => $file,
        output_file        => $newfile,
        tier_file_location => $tier_file_location,
    );
    unless ($tier_cmd->execute) {
        confess $self->error_message("Failed to tier variants successfully.");
    }
    return $newfile;
}

sub read_counts {
    my ($self, $variants_file, $directory) = @_;

    my $output_file = $self->result_file_path(
        input_file_path => $variants_file,
        suffix => 'rcnt',
        directory => $directory,
    );

    if (-s $variants_file) {
        my $bamfiles = join(",", $self->normal_bam, $self->tumor_bam);
        my $header = "Normal,Tumor";

        my $rc_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
            bam_files        => $bamfiles,
            output_file      => $output_file,
            variant_file     => $variants_file,
            genome_build     => $self->ref_seq_fasta,
            header_prefixes  => $header,
            indel_size_limit => 4,
        );
        unless ($rc_cmd->execute) {
            confess $self->error_message("Failed to obtain read counts for file $variants_file.");
        }
    }
    else {
        Genome::Sys->shellcmd( cmd => "touch $output_file" );
    }
    return $output_file;
}

sub annotation_build {
    my $self = shift;

    return $self->somatic_variation_build->annotation_build;
}

sub tumor_bam {
    my $self = shift;

    return $self->somatic_variation_build->tumor_build->merged_alignment_result->bam_path;
}

sub normal_bam {
    my $self = shift;

    return $self->somatic_variation_build->normal_build->merged_alignment_result->bam_path;
}

sub ref_seq_fasta {
    my $self = shift;

    return $self->somatic_variation_build->tumor_build->reference_sequence_build->full_consensus_path('fa');
}

sub snvs_dir {
    my $self = shift;
    return get_or_create_directory(
        File::Spec->join($self->_output_dir, 'snvs'),
    );
}

sub indels_dir {
    my $self = shift;
    return get_or_create_directory(
        File::Spec->join($self->_output_dir, 'indels'),
    );
}

sub review_dir {
    my $self = shift;
    return get_or_create_directory(
        File::Spec->join($self->_output_dir, 'review'),
    );
}

sub get_or_create_directory {
    my $path = shift;

    unless (-d $path) {
        Genome::Sys->create_directory($path);
    }

    return $path;
}

sub stage_snv_file {
    my $self = shift;

    $self->ensure_snvs_were_processed();

    my $snv_file = File::Spec->join($self->snvs_dir, 'snvs.hq.bed');
    my $cmd = sprintf('cat %s %s | joinx sort -o %s',
        File::Spec->join($self->_build_dir, 'effects', 'snvs.hq.novel.tier*.v2.bed'),
        File::Spec->join($self->_build_dir, 'effects', 'snvs.hq.previously_detected.tier*.v2.bed'),
        $snv_file,
    );
    Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => [$snv_file],
    );

    return $snv_file;
}

sub ensure_snvs_were_processed {
    my $self = shift;

    my $snv_file = File::Spec->join($self->_build_dir, 'effects', 'snvs.hq.novel.tier1.v2.bed');
    unless (-e $snv_file) {
        confess $self->error_message("SNV results for %s not found at %s", $self->sample_name, $snv_file);
    }
}

sub stage_indel_file {
    my $self = shift;

    $self->ensure_indels_were_processed();

    my $indel_file = File::Spec->join($self->indels_dir, 'indels.hq.bed');
    my $cmd = sprintf('cat %s %s | joinx sort -o %s',
        File::Spec->join($self->_build_dir, 'effects', 'indels.hq.novel.tier*.v2.bed'),
        File::Spec->join($self->_build_dir, 'effects', 'indels.hq.previously_detected.tier*.v2.bed'),
        $indel_file,
    );
    Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => [$indel_file],
    );

    return $indel_file;
}

sub ensure_indels_were_processed {
    my $self = shift;

    my $indel_file = File::Spec->join($self->_build_dir, 'effects', 'indels.hq.novel.tier1.v2.bed');
    unless (-e $indel_file) {
        confess $self->error_message("INDEL results for %s not found at %s", $self->sample_name, $indel_file);
    }
}

sub _filter_off_target_regions {
    my $self         = shift;
    my $file         = shift;
    my $directory    = shift;

    $self->status_message("Filtering out off-target regions for $file");

    my $featurelist = $self->get_or_create_featurelist_file();
    if (defined($featurelist)) {
        my $new_file = $self->result_file_path(
            input_file_path => $file,
            suffix => 'ontarget',
            directory => $directory,
        );

        Genome::Sys->shellcmd(cmd => "joinx intersect -a $file -b $featurelist -o $new_file");
        return $new_file;
    }
    else {
        $self->warning_message("feature list not found or target regions not specified; No target region filtering being done even though --restrict-to-target-regions set.");
    }

    return $file;
}

sub get_or_create_featurelist_file {
    my $self = shift;

    my $build = $self->somatic_variation_build;

    my $featurelist_file = File::Spec->join($self->_output_dir, 'featurelist');
    if (-e $featurelist_file) {
        return $featurelist_file;
    }

    my $featurelist;
    if ($self->target_regions) {
       $featurelist = $self->target_regions;
    }
    elsif ($build->tumor_build->can('target_region_set_name') and defined($build->tumor_build->target_region_set_name)) {
        my $featurelist_name = $build->tumor_build->target_region_set_name;
        $featurelist = Genome::FeatureList->get(name=>$featurelist_name)->file_path;
    }
    if (defined($featurelist) && (-s $featurelist)) {
        #clean up feature list
        my $outFh = Genome::Sys->open_file_for_writing("$featurelist_file.tmp");
        my $inFh  = Genome::Sys->open_file_for_reading($featurelist);
        while (my $line = $inFh->getline) {
            chomp($line);
            next if $line =~ /^track/;
            my ( $chr, $start, $stop, @rest) = split( /\t/, $line );
            #remove chr if present
            $chr =~ s/^chr//g;
            print $outFh join("\t",( $chr, $start, $stop, @rest)) . "\n";
        }
        close($inFh);
        close($outFh);

        Genome::Sys->shellcmd(cmd => "joinx sort $featurelist_file.tmp >$featurelist_file");
        Genome::Sys->shellcmd(cmd => "rm -f $featurelist_file.tmp");
        return $featurelist_file;
    }
    else {
        return;
    }
}

sub annotated_snvs_vcf {
    my $self = shift;

    return File::Spec->join($self->_build_dir, 'variants', 'snvs.annotated.vcf.gz');
}

sub _add_dbsnp_and_gmaf {
    my ($self, $snv_file, $indel_file) = @_;

    if (-s $self->annotated_snvs_vcf) {
        my $new_snv_file = $self->_add_dbsnp_and_gmaf_to_snv($snv_file);
        my $new_indel_file = $self->_add_dbsnp_and_gmaf_to_indel($indel_file);
        return ($new_snv_file, $new_indel_file);
    }
    else {
        $self->warning_message("Warning: couldn't find annotated SNV file in build, skipping dbsnp anno");
        return ($snv_file, $indel_file);
    }
}

sub _add_dbsnp_and_gmaf_to_snv {
    my $self       = shift;
    my $snv_file   = shift;

    my $output_file = $self->result_file_path(
        input_file_path => $snv_file,
        suffix => 'rsid',
        directory => $self->snvs_dir,
    );

    my $db_cmd = Genome::Model::Tools::Annotate::AddRsid->create(
        anno_file   => $snv_file,
        output_file => $output_file,
        vcf_file    => $self->annotated_snvs_vcf,
    );
    unless ($db_cmd->execute) {
        confess $self->error_message("Failed to add dbsnp anno to file $snv_file.");
    }
    return $output_file;
}

# pad indel file with tabs to match - if we ever start annotating with indels from dbsnp, replace this section
sub _add_dbsnp_and_gmaf_to_indel {
    my $self       = shift;
    my $indel_file = shift;

    my $output_file = $self->result_file_path(
        input_file_path => $indel_file,
        suffix => 'rsid',
        directory => $self->indels_dir,
    );

    my $outFh = Genome::Sys->open_file_for_writing($output_file);
    my $inFh  = Genome::Sys->open_file_for_reading($indel_file);
    while ( my $line = $inFh->getline ) {
        chomp($line);
        print $outFh $line . "\t\t\n"
    }
    close($inFh);
    close($outFh);

    return $output_file;
}

sub _create_master_files {
    my $self       = shift;
    my $snv_file   = shift;
    my $indel_file = shift;

    Genome::Sys->shellcmd(cmd => sprintf('tail -n +2 %s >> %s', $indel_file, $self->report . '.unsorted'));
    Genome::Sys->shellcmd(cmd => sprintf('tail -n +2 %s >> %s', $snv_file, $self->report . '.unsorted'));
    Genome::Sys->shellcmd(cmd => sprintf('joinx sort -i %s -o %s', $self->report . '.unsorted', $self->report . '.sorted'));

    # have to put the header on after joinx sort because joinx
    # won't recognize it as a bed format with a header
    Genome::Sys->shellcmd(cmd => sprintf('cat <(head -n 1 %s) %s > %s', $snv_file, $self->report . '.sorted', $self->report));
    Genome::Sys->shellcmd(cmd => sprintf('rm -f %s %s', $self->report . '.unsorted', $self->report . '.sorted'));

    # convert master table to excel
    my $workbook  = Spreadsheet::WriteExcel->new($self->report_xls);
    my $worksheet = $workbook->add_worksheet();

    my $row=0;
    my $inFh = Genome::Sys->open_file_for_reading($self->report);
    while (my $line = $inFh->getline) {
        chomp($line);
        my @F = split("\t", $line);
        for( my $i=0;$ i<@F; $i++) {
            $worksheet->write($row, $i, $F[$i]);
        }
        $row++;
    }
    close($inFh);
    $workbook->close();

    return 1;
}

sub _create_review_files {
    my $self = shift;

    $self->status_message("Generating Review files");

    $self->_create_review_bed();

    $self->_create_review_xml();
    $self->status_message("\n--------------------------------------------------------------------------------");
    $self->status_message("Sites to review are here:");
    $self->status_message($self->review_bed);
    $self->status_message("IGV XML file is here:");
    $self->status_message($self->review_xml);

    return 1;
}

sub _create_review_bed {
    my $self = shift;

    my @tiers = split(",", $self->tiers_to_review);
    my $tier_restricted_file = sprintf('%s.tier%s',
        $self->report, join('', @tiers));

    for my $i (@tiers) {
        Genome::Sys->shellcmd(
            cmd => sprintf('grep -w tier%s %s >> %s',
                $i, $self->report, $tier_restricted_file . '.tmp'),
        );
    }
    Genome::Sys->shellcmd(
        cmd => sprintf('joinx sort -i %s -o %s',
            $tier_restricted_file . '.tmp', $tier_restricted_file),
    );
    Genome::Sys->shellcmd(cmd => sprintf('rm -f %s', $tier_restricted_file . '.tmp'));
    annoFileToSlashedBedFile($tier_restricted_file, $self->review_bed);
}

sub _create_review_xml {
    my $self = shift;

    my $bam_files = join(",",($self->normal_bam, $self->tumor_bam));

    my $labels = sprintf('normal %s,tumor %s',
        $self->sample_name, $self->sample_name);

    #create the xml file for review
    my $dumpXML = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
        bams            => $bam_files,
        labels          => $labels,
        output_file     => $self->review_xml,
        genome_name     => $self->sample_name,
        review_bed_file => $self->review_bed,
        reference_name  => $self->igv_reference_name,
    );
    unless ($dumpXML->execute) {
        confess $self->error_message("Failed to create IGV xml file");
    }
}

1;
