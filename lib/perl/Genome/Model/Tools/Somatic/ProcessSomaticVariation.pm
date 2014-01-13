package Genome::Model::Tools::Somatic::ProcessSomaticVariation;

use warnings;
use strict;
use IO::File;
use Genome;
use Sort::Naturally qw(nsort);
use Genome::Info::IUB;
use Spreadsheet::WriteExcel;
use File::Slurp qw(read_dir);
use File::Basename;
use Carp qw(confess);

class Genome::Model::Tools::Somatic::ProcessSomaticVariation {
    is => 'Command::V2',
    has => [
        full_output_dir => {
            is => 'Text',
            is_optional => 1,
            is_transient => 1,
            is_mutable => 0,
            calculate => q { return $self->output_dir . "/" . $self->sample_name_dir; }
        },
        sample_name_dir => {
            is => 'Text',
            is_optional => 1,
            is_transient => 1,
            is_mutable => 0,
            calculate => q { return $self->get_unique_sample_name_dir(); }
        },
    ],
    has_input => [
        somatic_variation_model => {
            is => 'Genome::Model::SomaticVariation',
            doc => "SomaticVariation model",
        },
        output_dir => {
            is => 'Text',
            doc => "Directory where output will be stored (under a subdirectory with the sample name)",
        },
    ],
    has_optional_input => [
        igv_reference_name =>{
            is => 'Text',
            is_optional => 1,
            doc => "name of the igv reference to use",
            example_values => ["reference_build36","b37","mm9"],
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
        get_readcounts =>{
            is => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => "add readcounts to the final variant list",
        },
        restrict_to_target_regions =>{
            is => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => "only keep snv calls within the target regions. These are pulled from the build if possible",
        },
        target_regions =>{
            is => 'String',
            is_optional => 1,
            doc => "path to a target file region. Used in conjunction with --restrict-to-target-regions to limit sites to those appearing in these regions",
        },
        add_tiers =>{
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => "add tier information to the output",
        },
        variant_list_is_one_based => {
            is => 'Boolean',
            is_optional => 1,
            doc => "The variant list you provide is in annotation (one-based) format, instead of bed. This flag fixes that.",
            default => 0,
        },
        add_dbsnp_and_gmaf => {
            is => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => "if this is a recent build with vcf files (Jan 2013 or later), will add the rsids and GMAF information for all SNVs",
        },
        process_svs => {
            is => 'Boolean',
            is_optional => 1,
            doc => "Build has sv calls (probably WGS data). Most exomes won't have this",
            default => 0,
        },
        sites_to_pass => {
            is => 'String',
            is_optional => 1,
            doc => "an annotation file (5 col, 1 based) containing sites that will automatically be passed. This is useful when sequencing a relapse - the sites already found in the tumor don't need to be manually reviewed",
        },
        create_review_files => {
            is => 'Boolean',
            is_optional => 1,
            doc => "create xml and bed files for manual review",
            default => 0,
        },
        create_archive => {
            is => 'Boolean',
            is_optional => 1,
            doc => "create an archive suitable for passing to collaborators",
            default => 0,
        },
        include_vcfs_in_archive => {
            is => 'Boolean',
            is_optional => 1,
            doc => "include full vcf files in archive (very large files)",
            default => 0,
        },
        required_snv_callers => {
            is => 'Number',
            is_optional => 1,
            doc => "Number of independent algorithms that must call a SNV. If set to 1 (default), all calls are used",
            default => 1,
        },
        tiers_to_review => {
            is => 'String',
            is_optional => 1,
            doc => "comma-separated list of tiers to include in review. (e.g. 1,2,3 will create bed files with tier1, tier2, and tier3)",
            default => 1,
        },
        sample_name =>{
            is => 'Text',
            is_optional => 1,
            is_mutable => 0,
            calculate => q { return $self->somatic_variation_model->subject->name; },
            doc => "override the sample name on the build and use this name instead",
        }
    ],
};


sub help_detail {
  return <<HELP;
Given a SomaticVariation model, this tool will gather the resulting variants, remove
off-target sites, tier the variants, optionally filter them, etc. Calls are prepped for
manual review in the review/ directory.
HELP
}

sub _doc_authors {
  return <<AUTHS;
 Chris Miller
AUTHS
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
        $stop = $stop++;
    }
    else { #indel DEL or SNV
        $start = $start++;
    }
    return join("\t", ($chr, $start, $stop, $ref, $var));
}


sub annoToBed{
    my ($chr, $start, $stop, $ref, $var) = split("\t", $_[0]);
    if ($ref =~ /^[-*0]/) { #indel INS
        $stop = $stop--;
    }
    else { #indel DEL or SNV
        $start = $start--;
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
sub cleanFile{
    my ($file) = @_;

    my $newfile = addName($file,"clean");

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


sub getFilterSites{
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


sub removeFilterSites{
    my ($file, $filterSites) = @_;

    my $newfile = addName($file,"filtered");
    #handle zero size files
    if (-z $file) {
        Genome::Sys->shellcmd( cmd => "touch $newfile" );
        return $newfile;
    }

    my $outFh = Genome::Sys->open_file_for_writing($newfile);
    my $inFh  = Genome::Sys->open_file_for_reading($file);
    while (my $line = $inFh->getline) {
        chomp($line);
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
        if ($ref =~ /\//) {
            ( $ref, $var ) = split(/\//, $ref);
        }
        unless (defined($filterSites->{join("\t",($chr, $start, $stop, $ref, $var ))})) {
            print $outFh $line . "\n";
        }
    }
    close($outFh);
    return $newfile;
}



# sub getDeepestSubDir{ #danger - assumes one subdirectory per folder
#     my $dir = shift;
#     print STDERR "$dir\n";
#     if(grep { -d "$dir/$_"} read_dir($dir)){ #subdirectories exist
#         #get subdir, recurse into it
#         my @z = grep { -d "$dir/$_"} read_dir($dir);
#         @z = grep{ ! /indel/ } @z;
#         return(getDeepestSubDir("$dir/$z[0]"));
#     } else {
#         return($dir);
#     }
# }


##TODO - these two functions are brittle. There needs to be a better way to grab calls for specific callers. Ideally, from the VCF...
sub getVcfFile{
    my $self      = shift;
    my $build_dir = shift;
    my $dir       = shift;

    my $prefix = $build_dir . "/variants/snv/";

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


sub removeUnsupportedSites{
    my ($self, $snv_file) = @_;

    my $numcallers = $self->required_snv_callers;
    my $build_dir  = $self->somatic_variation_model->last_succeeded_build->data_directory;

    #hash all of the sites
    my $sites = $self->getFilterSites($snv_file);
    $self->status_message("here");
    for my $k (keys(%{$sites})) {
        $sites->{$k} = 0;
    }


    #Look for the callers
    my @dirs = map {basename($_) } glob("$build_dir/variants/snv/*");
    #remove non-caller dirs
    @dirs = grep{ ! /^intersect|^union|^samtools/ } @dirs;

    #count the number of callers that called each site from the vcfs
    for my $dir (@dirs) {
        my $file = $self->getVcfFile($build_dir, $dir);
        my $ifh = Genome::Sys->open_gzip_file_for_reading($file);

        while (my $line = $ifh->getline) {
            chomp $line;
            next if $line =~ /^#/;
            my ($chr, $pos, $id, $ref, $var, @rest) = split /\t/, $line;
            my @vars = split(",", $var);
            for my $v (@vars) {
                my $key = join("\t",($chr, $pos--, $pos, $ref, $v));
                if (defined($sites->{$key})) {
                    $sites->{$key} = $sites->{$key}++;
                }
            }
        }
        $ifh->close;
    }

    my $ofh = Genome::Sys->open_file_for_writing(addName($snv_file, "gt" . $numcallers . "callers"));
    #read the snv_file again to preserve order, traiing fields, etc.
    my $ifh = Genome::Sys->open_file_for_reading($snv_file);
    while (my $line = $ifh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $ref, $var, @rest) = split /\t/, $line;
        my $key = join("\t",($chr, $start, $stop, $ref, $var));

        if (!defined($sites->{$key})) {
            $self->status_message("wut?: " . $key);
        }
        if ($sites->{$key} >= $numcallers) {
            print $ofh $line . "\n";
        }
    }
    close($ifh);
    close($ofh);
    return addName($snv_file, "gt" . $numcallers . "callers");
}




sub doAnnotation{
    my $self                  = shift;
    my $file                  = shift;
    my $outfile               = shift;

    my $annotation_build_name = $self->somatic_variation_model->annotation_build->name;

    if ($file =~ /.bed$/) {
        $file = bedFileToAnnoFile($file);
    }

    my $newfile;
    if ($outfile) {
        $newfile = $outfile;
    }
    else {
        $newfile = $file . ".anno";
    }

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


sub addTiering{
    my $self               = shift;
    my $file               = shift;
    my $tier_file_location = shift;
    my $outfile            = shift;

    my $newfile = addName($file, "tiered");

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

sub getReadcounts{
    my ($self, $file, $ref_seq_fasta, @bams) = @_;
    #todo - should check if input is bed and do coversion if necessary

    if (-s "$file") {
        my $bamfiles = join(",",@bams);
        my $header = "Tumor";
        if (@bams == 2) {
            $header = "Normal,Tumor";
        }
        #get readcounts from the tumor bam only
        my $rc_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
            bam_files        => $bamfiles,
            output_file      => "$file.rcnt",
            variant_file     => "$file",
            genome_build     => $ref_seq_fasta,
            header_prefixes  => $header,
            indel_size_limit => 4,
        );
        unless ($rc_cmd->execute) {
            confess $self->error_message("Failed to obtain readcounts for file $file.");
        }
    }
    else {
        Genome::Sys->shellcmd( cmd => "touch $file.rcnt" );
    }
    return "$file.rcnt";
}

sub get_unique_sample_name_dir {
    my $self = shift;

    my $output_dir  = $self->output_dir;
    my $sample_name = $self->sample_name;

    if (-e "$output_dir/$sample_name") {
        my $suffix = 1;
        my $newname = $sample_name . "-" . $suffix;
        while (-e "$output_dir/$newname") {
            $suffix++;
            $newname = $sample_name . "-" . $suffix;
        }
        return $newname;
    }
    else {
        return $sample_name;
    }
}

sub create_directories {
    my $self = shift;

    my $full_output_dir = $self->full_output_dir;
    my $output_dir      = $self->output_dir;
    my $build_dir       = $self->somatic_variation_model->last_succeeded_build->data_directory;

    Genome::Sys->create_directory("$full_output_dir");
    Genome::Sys->create_directory("$full_output_dir/snvs");
    Genome::Sys->create_directory("$full_output_dir/indels");
    unless ($self->create_review_files) {
        Genome::Sys->create_directory("$output_dir/review");
    }
    Genome::Sys->create_symlink($build_dir, "$full_output_dir/build_directory");

    return 1;
}

sub stage_snv_file {
    my $self = shift;

    my $full_output_dir = $self->full_output_dir;
    my $build_dir       = $self->somatic_variation_model->last_succeeded_build->data_directory;

    my $snv_file = "$build_dir/effects/snvs.hq.novel.tier1.v2.bed";
    unless (-e $snv_file) {
        confess $self->error_message("SNV results for " . $self->sample_name . " not found at $snv_file");
    }

    #cat all the filtered snvs together
    $snv_file = "$full_output_dir/snvs/snvs.hq.bed";
    Genome::Sys->shellcmd( cmd => "cat $build_dir/effects/snvs.hq.novel.tier*.v2.bed $build_dir/effects/snvs.hq.previously_detected.tier*.v2.bed | joinx sort >$snv_file");

#  `ln -s $snv_file $full_output_dir/snvs/` unless( -e "$full_output_dir/snvs/$snv_file");

    return $snv_file;
}

sub stage_indel_file {
    my $self = shift;

    my $full_output_dir = $self->full_output_dir;
    my $build_dir       = $self->somatic_variation_model->last_succeeded_build->data_directory;

    my $indel_file = "$build_dir/effects/indels.hq.novel.tier1.v2.bed";
    unless (-e $indel_file) {
        confess $self->error_message("INDEL results for " . $self->sample_name . " not found at $indel_file");
    }

    #cat all the filtered indels together
    $indel_file = "$full_output_dir/indels/indels.hq.bed";
    Genome::Sys->shellcmd( cmd => "cat $build_dir/effects/indels.hq.novel.tier*.v2.bed $build_dir/effects/indels.hq.previously_detected.tier*.v2.bed | joinx sort >$indel_file" );

#  `ln -s $indel_file $full_output_dirindels/` unless( -e "$full_output_dir/indels/$indel_file");

    return $indel_file;
}

sub stage_sv_file {
    my $self = shift;

    my $full_output_dir = $self->full_output_dir;
    my $build_dir       = $self->somatic_variation_model->last_succeeded_build->data_directory;

    my $sv_file;
#THIS NEEDS A TEST TO MAKE SURE THAT $self->process_svs IS BEING SET CORRECTLY
    if ($self->process_svs) {
        my @sv_files = glob("$build_dir/variants/sv/union-union-sv_breakdancer_*sv_squaredancer*/svs.merge.file.somatic");
        $sv_file = $sv_files[0];
        unless (-e $sv_file) {
            $self->warning_message("SV results for " . $self->sample_name . " not found, skipping SVs");
            $self->process_svs(0);
        }

        Genome::Sys->create_directory("$full_output_dir/svs");
        Genome::Sys->create_symlink($sv_file, "$full_output_dir/svs/svs.hq") unless( -e "$full_output_dir/svs/$sv_file");

#       #annotate the svs
#       my $anno_cmd = Genome::Model::Tools::Annotate::Sv::Combine->create(
#           input-file => $sv_file,
#           output-file => "$full_output_dir/svs/svs.hq.annotated",
#           annotation-build
#           dbsnp-annotation-file /gsc/scripts/opt/genome/db/genome-db-dbsnp/human/build/37/132/dbsnp.csv
#           dbvar-annotation-file /gsc/scripts/opt/genome/db/dbvar/human/build37/dbvar.tsv
#           fusion-transcripts-fusion-output-file /tmp/zout.fusions
#           repeat-masker-annotation-file /gscuser/aregier/git/genome/vep/lib/perl/Genome/repeat_masker.tsv
#           annotator-list=Transcripts,FusionTranscripts,Dbsnp,Segdup,Dbvar
#           segdup-annotation-file /gsc/scripts/opt/genome/db/ucsc/human/build37/segdup.tsv
#           chrA-column 1
#           bpA-column 2
#           chrB-column 4
#           bpB-column 5
#           event-type-column 7
#           score-column 12
#           orient-column 8

#           );
#       unless ($anno_cmd->execute) {
#           die "Failed to annotate sv file\n";
#       }
# #gmt annotate sv combine --input-file /tmp/zin --output-file /tmp/zout1 --annotation-build 124434505 --dbsnp-annotation-file /gsc/scripts/opt/genome/db/genome-db-dbsnp/human/build37/132/dbsnp.csv --dbvar-annotation-file /gsc/scripts/opt/genome/db/dbvar/human/build37/dbvar.tsv --fusion-transcripts-fusion-output-file /tmp/zout.fusions --repeat-masker-annotation-file /gscuser/aregier/git/genome/vep/lib/perl/Genome/repeat_masker.tsv --annotator-list=Transcripts,FusionTranscripts,Dbsnp,Segdup,Dbvar --segdup-annotation-file /gsc/scripts/opt/genome/db/ucsc/human/build37/segdup.tsv --chrA-column 1 --bpA-column 2 --chrB-column 4 --bpB-column 5 --event-type-column 7 --score-column 12 --orient-column 8
    }

    return $sv_file;
}

#########################################################################################################
sub execute {
    my $self = shift;

    my $output_dir = $self->output_dir;
    $output_dir =~ s/(\/)+$//; # Remove trailing forward-slashes if any

    unless (-e $output_dir) {
        Genome::Sys->create_directory($output_dir);
    }

    # Check on the input data before starting work
    my $model = $self->somatic_variation_model;

    unless (-e $output_dir) {
        confess $self->error_message("Output directory not found: $output_dir");
    }

    #grab the info from the model
    my $build = $model->last_succeeded_build;
    unless (defined($build)) {
        confess $self->error_message("Model " . $model->id . "has no succeeded builds");
    }

    my $tumor_model = $model->tumor_model;
    my $normal_model = $model->normal_model;

    my $ref_seq_build = $tumor_model->reference_sequence_build;
    my $ref_seq_fasta = $ref_seq_build->full_consensus_path('fa');
    my $tiering_files = $model->annotation_build->data_directory . "/annotation_data/tiering_bed_files_v3/";
    $self->status_message("Processing model with sample_name: " . $self->sample_name);

    my $tumor_bam = $tumor_model->last_succeeded_build->merged_alignment_result->bam_file;
    my $normal_bam = $normal_model->last_succeeded_build->merged_alignment_result->bam_file;
    my $build_dir = $build->data_directory;

    my $igv_reference_name = $self->igv_reference_name;
    if ($self->create_review_files && !defined($self->igv_reference_name)) {
        confess $self->error_message("igv-reference-name required if --create-review-files is specified");
    }


    # create subdirectories, get files in place

    # if multiple models with the same name, add a suffix
    my $sample_name_dir = $self->sample_name_dir;
    my $full_output_dir = $self->full_output_dir;

    #make the directory structure
    $self->create_directories();

    # Check if the necessary files exist in this build and put them in the processing location
    my $snv_file   = $self->stage_snv_file();
    my $indel_file = $self->stage_indel_file();
    my $sv_file    = $self->stage_sv_file();

    #--------------------------------------------------------------
    # munge through SNV file to remove duplicates and fix IUB codes
    $snv_file = cleanFile($snv_file);
    $indel_file = cleanFile($indel_file);

    #-------------------------------------------------
    # filter out the off-target regions, if target regions are available
    if ($self->restrict_to_target_regions) {
        $snv_file = $self->_filter_off_target_regions($snv_file);
        $indel_file = $self->_filter_off_target_regions($indel_file);
    }

    #-------------------------------------------------
    # remove filter sites specified by the user
    if (defined($self->filter_sites)) {
        $self->status_message("Applying user-supplied filter");
        my $filterSites = $self->getFilterSites($self->filter_sites);
        $snv_file   = removeFilterSites($snv_file, $filterSites);
        $indel_file = removeFilterSites($indel_file, $filterSites);
    }

    #-------------------------------------------------
    # remove filter regions specified by the user
    if (defined($self->filter_regions)) {
        $snv_file = $self->_filter_regions($snv_file);
        $indel_file = $self->_filter_regions($indel_file);
    }

    #-------------------------------------------------------
    # remove regions called by less than the required number of callers
    unless ($self->required_snv_callers == 1) {
        $self->status_message("Removing regions supported by less than " . $self->required_snv_callers . " regions");
        $snv_file = $self->removeUnsupportedSites($snv_file);
    }

    #------------------------------------------------------
    # do annotation
    $snv_file   = $self->doAnnotation($snv_file);
    $indel_file = $self->doAnnotation($indel_file);

    #-------------------------------------------------------
    # add tiers
    if ($self->add_tiers) {
        $self->status_message("Adding tiers");
        #do annotation
        $snv_file   = $self->addTiering($snv_file, $tiering_files);
        $indel_file = $self->addTiering($indel_file, $tiering_files);

        #convert back to annotation format (1-based)
        $snv_file = bedFileToAnnoFile($snv_file);
        $indel_file = bedFileToAnnoFile($indel_file);
    }

    #----------------------------------------------------
    # add dbsnp/gmaf
    if ($self->add_dbsnp_and_gmaf) {
        $snv_file   = $self->_add_dbsnp_and_gmaf_to_snv($build_dir, $snv_file);
        $indel_file = $self->_add_dbsnp_and_gmaf_to_indel($build_dir, $indel_file);
    }

    #-------------------------------------------------------
    # get readcounts
    if ($self->get_readcounts) {
      $self->status_message("Getting readcounts");
      if (-s "$snv_file") {
          $snv_file   = $self->getReadcounts($snv_file, $ref_seq_fasta, $normal_bam, $tumor_bam);
      }
      if (-s "$indel_file") {
          $indel_file = $self->getReadcounts($indel_file, $ref_seq_fasta, $normal_bam, $tumor_bam);
      }
    }

    #------------------------------------------------------
    # combine the files into one master table
    $self->_create_master_files($snv_file, $indel_file);

    #------------------------------------------------------
    # now get the files together for review
    if ($self->create_review_files) {
        $self->_create_review_files();
    }

    #------------------------------------------------
    # tar up the files to be sent to collaborators
    if ($self->create_archive) {
        $self->_create_archive($build_dir, $sv_file);
    }

    return 1;
}

sub _filter_off_target_regions {
    my $self = shift;
    my $file = shift;

    $self->status_message("Filtering out off-target regions for $file");

    my $featurelist = $self->get_or_create_featurelist_file();
    if (defined($featurelist)) {
          my $new_file = addName($file,"ontarget");

          Genome::Sys->shellcmd(cmd => "joinx intersect -a $file -b $featurelist >$new_file");
          $file = $new_file;
    }
    else {
        $self->warning_message("feature list not found or target regions not specified; No target region filtering being done even though --restrict-to-target-regions set.");
    }

    return $file;
}

sub get_or_create_featurelist_file {
    my $self = shift;

    my $full_output_dir = $self->full_output_dir;
    my $model = $self->somatic_variation_model;

    my $featurelist_file = "$full_output_dir/featurelist";

    if (-e $featurelist_file) {
        return $featurelist_file;
    }

    my $featurelist;
    if ($self->target_regions) {
       $featurelist = $self->target_regions;
    }
    elsif ($model->tumor_model->can('target_region_set_name') and defined($model->tumor_model->target_region_set_name)) {
        my $featurelist_name = $model->tumor_model->target_region_set_name;
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
        return 0;
    }
}

sub _filter_regions {
    my $self = shift;
    my $file = shift;

    my $filter_file = $self->get_or_create_filter_file();

    $self->status_message("Removing user-specified filter for $file");
    my $filtered_file = "$file.filteredReg";
    Genome::Sys->shellcmd( cmd => "joinx intersect --miss-a $filtered_file -a $file -b $filter_file >/dev/null" );

    return $filtered_file;
}

sub get_or_create_filter_file {
    my $self = shift;

    my $full_output_dir = $self->full_output_dir;

    my @filters = split(",", $self->filter_regions);

    my $filter_file = "$full_output_dir/filter";

    if (-e $filter_file) {
        return $filter_file;
    }

    Genome::Sys->concatenate_files(\@filters, $filter_file);

    return $filter_file;
}

sub _add_dbsnp_and_gmaf_to_snv {
    my $self       = shift;
    my $build_dir  = shift;
    my $snv_file   = shift;

    $self->status_message("==== adding dbsnp ids ====");
    $self->status_message("$build_dir/variants/snvs.annotated.vcf.gz");
    if (-s "$build_dir/variants/snvs.annotated.vcf.gz") {
        my $db_cmd = Genome::Model::Tools::Annotate::AddRsid->create(
            anno_file   => $snv_file,
            output_file => "$snv_file.rsid",
            vcf_file    => "$build_dir/variants/snvs.annotated.vcf.gz",
        );
        unless ($db_cmd->execute) {
            confess $self->error_message("Failed to add dbsnp anno to file $snv_file.");
        }
        $snv_file = "$snv_file.rsid";
        return $snv_file;
    }
    else {
        $self->warning_message("Warning: couldn't find annotated SNV file in build, skipping dbsnp anno");
    }
}

sub _add_dbsnp_and_gmaf_to_indel {
    my $self       = shift;
    my $build_dir  = shift;
    my $indel_file = shift;

    $self->status_message("==== padding indel file with tabs to match ====");
    $self->status_message("$build_dir/variants/snvs.annotated.vcf.gz");
    if (-s "$build_dir/variants/snvs.annotated.vcf.gz") {
        #pad indel file with tabs to match - if we ever start annotating with indels from dbsnp, replace this section
        my $outFh = Genome::Sys->open_file_for_writing("$indel_file.rsid");
        my $inFh  = Genome::Sys->open_file_for_reading($indel_file);
        while ( my $line = $inFh->getline ) {
            chomp($line);
            print $outFh $line . "\t\t\n"
        }
        close($inFh);
        close($outFh);

        $indel_file = "$indel_file.rsid";
        return $indel_file;
    }
    else {
        $self->warning_message("Couldn't find annotated SNV file in build, skipping dbsnp anno");
    }
}

sub _create_master_files {
    my $self       = shift;
    my $snv_file   = shift;
    my $indel_file = shift;

    my $full_output_dir = $self->full_output_dir;

    Genome::Sys->shellcmd(cmd => "head -n 1 $snv_file >$full_output_dir/snvs.indels.annotated");
    Genome::Sys->shellcmd(cmd => "tail -n +2 $indel_file >>$full_output_dir/snvs.indels.annotated.tmp");
    Genome::Sys->shellcmd(cmd => "tail -n +2 $snv_file >>$full_output_dir/snvs.indels.annotated.tmp");
    Genome::Sys->shellcmd(cmd => "joinx sort -i $full_output_dir/snvs.indels.annotated.tmp >>$full_output_dir/snvs.indels.annotated");
    Genome::Sys->shellcmd(cmd => "rm -f $full_output_dir/snvs.indels.annotated.tmp");

    # convert master table to excel
    my $workbook  = Spreadsheet::WriteExcel->new("$full_output_dir/snvs.indels.annotated.xls");
    my $worksheet = $workbook->add_worksheet();

    my $row=0;
    my $inFh = Genome::Sys->open_file_for_reading("$full_output_dir/snvs.indels.annotated");
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

    my $review_dir      = $self->output_dir . "/review";
    my $sample_name_dir = $self->sample_name_dir;
    my $full_output_dir = $self->full_output_dir;

    $self->status_message("Generating Review files");
    my @tiers = split(",", $self->tiers_to_review);
    my $tierstring = join("",@tiers);
    for my $i (@tiers) {
        Genome::Sys->shellcmd(cmd => "grep -w tier$i $full_output_dir/snvs.indels.annotated >>$full_output_dir/snvs.indels.annotated.tier$tierstring.tmp");
    }
    Genome::Sys->shellcmd(cmd => "joinx sort -i $full_output_dir/snvs.indels.annotated.tier$tierstring.tmp >$full_output_dir/snvs.indels.annotated.tier$tierstring");
    Genome::Sys->shellcmd(cmd => "rm -f $full_output_dir/snvs.indels.annotated.tier$tierstring.tmp");
    annoFileToSlashedBedFile("$full_output_dir/snvs.indels.annotated.tier$tierstring","$review_dir/$sample_name_dir.bed");

    my $tumor_bam = $self->somatic_variation_model->tumor_model->last_succeeded_build->merged_alignment_result->bam_file;
    my $normal_bam = $self->somatic_variation_model->normal_model->last_succeeded_build->merged_alignment_result->bam_file;

    my $bam_files = join(",",($normal_bam, $tumor_bam));
    my $labels = join(",",("normal $sample_name_dir","tumor $sample_name_dir"));

    my $igv_reference_name = $self->igv_reference_name;
    unless (defined($igv_reference_name)) {
          $self->warning_message("No IGV reference name supplied - defaulting to build 37");
    }

    #create the xml file for review
    my $dumpXML = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
        bams            => "$bam_files",
        labels          => "$labels",
        output_file     => "$review_dir/$sample_name_dir.xml",
        genome_name     => $sample_name_dir,
        review_bed_file => "$review_dir/$sample_name_dir.bed",
        reference_name  => $igv_reference_name,
    );
    unless ($dumpXML->execute) {
        confess $self->error_message("Failed to create IGV xml file");
    }

    $self->status_message("\n--------------------------------------------------------------------------------");
    $self->status_message("Sites to review are here:");
    $self->status_message("$review_dir/$sample_name_dir.bed");
    $self->status_message("IGV XML file is here:");
    $self->status_message("$review_dir/$sample_name_dir.xml");

    return 1;
}

sub _create_archive {
    my $self        = shift;
    my $build_dir   = shift;
    my $sv_file     = shift;

    my $full_output_dir = $self->full_output_dir;
    my $archive_dir     = "$full_output_dir/" . $self->sample_name_dir;
    Genome::Sys->create_directory("$archive_dir");

    #symlink VCF files
    if ($self->include_vcfs_in_archive) {
        if (-e "$build_dir/variants/indels.detailed.vcf.gz") {
            Genome::Sys->create_symlink("$build_dir/variants/indels.detailed.vcf.gz", "$archive_dir/indels.vcf.gz");
        }
        elsif (-e "$build_dir/variants/indels.vcf.gz") {
            Genome::Sys->create_symlink("$build_dir/variants/indels.vcf.gz", "$archive_dir/indels.vcf.gz");
        }
        else {
            $self->warning_message("No indel VCF file available. If this is an older model, a rebuild may fix this");
        }

        if (-e "$build_dir/variants/snvs.annotated.vcf.gz") {
            Genome::Sys->create_symlink("$build_dir/variants/snvs.annotated.vcf.gz", "$archive_dir/snvs.vcf.gz");
        }
        elsif (-e "$build_dir/variants/snvs.vcf.gz") {
            Genome::Sys->create_symlink("$build_dir/variants/snvs.vcf.gz", "$archive_dir/snvs.vcf.gz");
        }
        else {
            $self->warning_message("No snv VCF file available. If this is an older model, a rebuild may fix this");
        }
    }

    #symlink annotated snvs and indels
    Genome::Sys->create_symlink("$full_output_dir/snvs.indels.annotated", "$archive_dir/snvsAndIndels.annotated");
    #symlink annoted snvs and indels excel file
    Genome::Sys->create_symlink("$full_output_dir/snvs.indels.annotated.xls", "$archive_dir/snvsAndIndels.annotated.xls");
    #symlink sv calls
    if ($self->process_svs) {
        Genome::Sys->create_symlink("$sv_file", "$archive_dir/svs");
    }
    #synlink cnv calls - todo

    #tar it up
#WRITE Genome::Sys TEST FOR THE NEW options PARAMETER
    Genome::Sys->tar(tar_path => "$archive_dir.tar.gz", input_directory => "$archive_dir", options => "-chzvf");

    return 1;
}

1;
