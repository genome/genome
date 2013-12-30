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

class Genome::Model::Tools::Somatic::ProcessSomaticVariation {
    is => 'Command',
    has => [
        somatic_variation_model => {
            is => 'Genome::Model::SomaticVariation',
            is_optional => 1,
            id_by => 'somatic_variation_model_id',
        },
    ],
    has_input => [
        somatic_variation_model_id => {
            is => 'Text',
            is_optional => 0,
            doc => "ID of SomaticVariation model",
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
            doc => "override the sample name on the build and use this name instead",
        },
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
    my ($file,$outfile) = @_;


    #remove bed from name
    my $newfile = $file;
    $newfile =~ s/\.bed//g;

    if($outfile){
        $newfile = $outfile;
    }

    open(OUTFILE,">$newfile");
    my $inFh = IO::File->new( $file ) || die "can't open file2\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        if($line =~ /^chrom/){
            print OUTFILE $line . "\n";
            next;
        }
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        if($ref =~ /\//){
            my @alleles = split("/",$ref);
            print OUTFILE bedToAnno(join("\t",($chr, $start, $stop, $alleles[0], $alleles[1]))) . "\t" . join("\t",($var,@rest)) . "\n";
        } else {
            print OUTFILE bedToAnno(join("\t",($chr, $start, $stop, $ref, $var))) . "\t" . join("\t",@rest) . "\n";
        }
    }
    close(OUTFILE);
    close($inFh);
    return($newfile);
}


sub annoFileToBedFile{
    my ($file,$outfile) = @_;
    #add bed to name
    my $newfile = $file . ".bed";
    if($outfile){
        $newfile = $outfile;
    }

    open(OUTFILE,">$newfile");
    my $inFh = IO::File->new( $file ) || die "can't open file2\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        if($line =~ /^chrom/){
            print OUTFILE $line . "\n";
            next;
        }
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        print OUTFILE annoToBed(join("\t",($chr, $start, $stop, $ref, $var))) . "\t" . join("\t",@rest) . "\n";
    }
    close(OUTFILE);
    close($inFh);
    return($newfile);
}

sub annoFileToSlashedBedFile{
    my ($file,$outfile) = @_;

    my $newfile = $file . ".bed";
    if($outfile){
        $newfile = $outfile;
    }

    open(OUTFILE,">$newfile");
    my $inFh = IO::File->new( $file ) || die "can't open file2\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        if($line =~ /^chrom/){
            next;
        }
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        my $bed = annoToBed(join("\t",($chr, $start, $stop, $ref, $var)));
        my @bedline = split(/\t/,$bed);
        print OUTFILE join("\t",(@bedline[0..2],"$bedline[3]/$bedline[4]")) . "\n";
    }
    close(OUTFILE);
    close($inFh);
    return($newfile);
}

sub bedToAnno{
    my ($chr,$start,$stop,$ref,$var) = split("\t",$_[0]);
    #print STDERR join("|",($chr,$start,$stop,$ref,$var)) . "\n";
    if ($ref =~ /^[-0*]/){ #indel INS
        $stop = $stop+1;
    } else { #indel DEL or SNV
        $start = $start+1;
    }
    return(join("\t",($chr,$start,$stop,$ref,$var)));
}


sub annoToBed{
    my ($chr,$start,$stop,$ref,$var) = split("\t",$_[0]);
    if ($ref =~ /^[-*0]/){ #indel INS
        $stop = $stop-1;
    } else { #indel DEL or SNV
        $start = $start-1;
    }
    #handle 5 col or 4 col ref/var
    if(defined($var)){
        return(join("\t",($chr,$start,$stop,$ref,$var)));
    } else {
        return(join("\t",($chr,$start,$stop,$ref)));
    }

}

sub intersects{
    my ($st,$sp,$st2,$sp2) = @_;
    if((($sp2 >= $st) && ($sp2 <= $sp)) ||
       (($sp >= $st2) && ($sp <= $sp2))){
        return 1;
    }
    return 0;
}

sub fixIUB{
    my ($ref,$var) = @_;
    my @vars = Genome::Info::IUB->variant_alleles_for_iub($ref,$var);
    return @vars;
}

sub addName{
    my ($file,$name) = @_;
    if($file =~ /\.bed$/){
        $file =~ s/\.bed//g;
        return($file . "." . $name . ".bed");
    } else {
        return($file . "." . $name);
    }
}

#read in the file, output a cleaned-up version
sub cleanFile{
    my ($file) = @_;

    my $newfile = addName($file,"clean");

    my %dups;
    open(OUTFILE, ">$newfile");

    my $inFh = IO::File->new( $file ) || die "can't open file $file\n";
    while( my $line = $inFh->getline ){
        chomp($line);
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        if($ref =~ /\//){
            ( $ref, $var ) = split(/\//, $ref);
        }

        $ref =~ s/0/-/g;
        $var =~ s/0/-/g;
        $ref =~ s/\*/-/g;
        $var =~ s/\*/-/g;

        my @vars = ($var);
        unless($ref =~ /-/ || $var =~ /-/){ #fixiub doesn't handle indels
            @vars = fixIUB($ref, $var);
        }

        foreach my $v (@vars){
            unless (exists($dups{join("\t",($chr, $start, $stop, $ref, $v ))})){
                print OUTFILE join("\t",($chr, $start, $stop, $ref, $v )) . "\n";
            }
            $dups{join("\t",($chr, $start, $stop, $ref, $v ))} = 1;
        }
    }
    close(OUTFILE);
    close($inFh);
    Genome::Sys->shellcmd( cmd => "joinx sort -i $newfile >$newfile.tmp" );
    Genome::Sys->shellcmd( cmd => "mv -f $newfile.tmp $newfile");
    return($newfile)
}


sub getFilterSites{
    my ($filter_string) = @_;
    my @filters = split(",",$filter_string);
    my %filterSites;

    foreach my $filterfile (@filters){
        if( -s $filterfile){
            #store sites to filter out in a hash
            my $inFh = IO::File->new( $filterfile ) || die "can't open file5\n";
            while( my $line = $inFh->getline )
            {
                chomp($line);
                #handle either 5 col (Ref\tVar) or 4 col (Ref/Var) bed
                my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
                if($ref =~ /\//){
                    ( $ref, $var ) = split(/\//, $ref);
                }
                $ref =~ s/0/-/g;
                $var =~ s/0/-/g;
                $ref =~ s/\*/-/g;
                $var =~ s/\*/-/g;

                my @vars = fixIUB($ref, $var);
                foreach my $v (@vars){
                    $filterSites{join("\t",($chr, $start, $stop, $ref, $v ))} = 0;
                }
            }
            close($inFh);
        } else {
            die("filter sites file doesn't exist or has zero size: " . $filterfile);
        }
    }
    return(\%filterSites);
}


sub removeFilterSites{
    my ($file,$filterSites) = @_;

    my $newfile = addName($file,"filtered");
    #handle zero size files
    if( -z $file ){
        Genome::Sys->shellcmd( cmd => "touch $newfile" );
        return($newfile);
    }

    open(FILFILE,">$newfile") || die ("couldn't open filter output file");
    my $inFh = IO::File->new( $file ) || die "can't open filter input file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
        if($ref =~ /\//){
            ( $ref, $var ) = split(/\//, $ref);
        }
        unless (defined($filterSites->join("\t",($chr, $start, $stop, $ref, $var )))){
            print FILFILE $line . "\n";
        }
    }
    close(FILFILE);
    return($newfile)
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
    my ($build_dir, $dir) = @_;
    my $prefix = $build_dir . "/variants/snv/";

    if($dir=~/strelka/){
        if(-s "$prefix/$dir/snvs.vcf.gz"){
            return("$prefix/$dir/snvs.vcf.gz");
        }

    } elsif($dir=~/varscan/){
        if(-s "$prefix/$dir/varscan-high-confidence-v1-d41d8cd98f00b204e9800998ecf8427e/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/snvs.vcf.gz"){
            return("$prefix/$dir/varscan-high-confidence-v1-d41d8cd98f00b204e9800998ecf8427e/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/snvs.vcf.gz");
        }

    } elsif($dir=~/sniper/){
        if(-s "$prefix/$dir/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/somatic-score-mapping-quality-v1-39b60f48b6f8c9e63436a5424305e9fd/snvs.vcf.gz"){
            return("$prefix/$dir/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/somatic-score-mapping-quality-v1-39b60f48b6f8c9e63436a5424305e9fd/snvs.vcf.gz");
        }

    } else {
        die("Don't know how to find the calls for ")
    }
    return("ERROR: couldn't find the snvs.vcf.gz file under directory $dir\n");
}


sub removeUnsupportedSites{
    my ($snv_file, $numcallers, $build_dir) = @_;
    
    #hash all of the sites
    my $sites = getFilterSites($snv_file);
    print "here\n";
    for my $k (keys(%{$sites})){
        $sites->{$k} = 0;
    }


    #Look for the callers
    my @dirs = map {basename($_) } glob("$build_dir/variants/snv/*");
    #remove non-caller dirs
    @dirs = grep{ ! /^intersect|^union|^samtools/ } @dirs;

    #count the number of callers that called each site from the vcfs
    for my $dir (@dirs){
        my $file = getVcfFile($build_dir,$dir);
        my $ifh = Genome::Sys->open_gzip_file_for_reading($file);

        while (my $line = $ifh->getline) {
            chomp $line;
            next if $line =~ /^#/;
            my ($chr, $pos, $id, $ref, $var, @rest) = split /\t/, $line;
            my @vars = split(",",$var);
            for my $v (@vars){
                my $key = join("\t",($chr, $pos-1, $pos, $ref, $v));
                if(defined($sites->{$key})){
                    $sites->{$key} = $sites->{$key}+1;
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
        
        if(!defined($sites->{$key})){
            print STDERR "wut?: " . $key . "\n";
        }
        if ($sites->{$key} >= $numcallers){
            print $ofh $line . "\n";
        }
    }
    close($ifh);
    close($ofh);
    return(addName($snv_file, "gt" . $numcallers . "callers"));
}




sub doAnnotation{
    my ($file,$annotation_build_name, $outfile) = @_;

    if($file =~ /.bed$/){
        $file = bedFileToAnnoFile($file);
    }

    my $newfile;
    if($outfile){
        $newfile = $outfile;
    }
    else {
        $newfile = $file . ".anno";
    }

    #handle zero size files
    if( -z $file ){
        Genome::Sys->shellcmd( cmd => "touch $newfile" );
        return($newfile);
    }

    my $anno_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
        variant_file => $file,
        output_file => $newfile,
        reference_transcripts => $annotation_build_name,
        annotation_filter => "top",
    );
    unless ($anno_cmd->execute) {
        die "Failed to annotate variants successfully.\n";
    }
    return($newfile);
}


sub addTiering{
    my ($file, $tier_file_location, $outfile) = @_;

    my $newfile = addName($file, "tiered");

    #handle zero size files
    if( -z $file ){
        Genome::Sys->shellcmd( cmd => "touch $newfile" );
        return($newfile);
    }

    my $tier_cmd = Genome::Model::Tools::FastTier::AddTiers->create(
        input_file => $file,
        output_file => $newfile,
        tier_file_location => $tier_file_location,
    );
    unless ($tier_cmd->execute) {
        die "Failed to tier variants successfully.\n";
    }
    return($newfile);
}

sub getReadcounts{
    my ($file, $ref_seq_fasta, @bams) = @_;
    #todo - should check if input is bed and do coversion if necessary

    if( -s "$file" ){
        my $bamfiles = join(",",@bams);
        my $header = "Tumor";
        if(@bams == 2){
            $header = "Normal,Tumor";
        }
        #get readcounts from the tumor bam only
        my $rc_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
            bam_files => $bamfiles,
            output_file => "$file.rcnt",
            variant_file => "$file",
            genome_build => $ref_seq_fasta,
            header_prefixes => $header,
            indel_size_limit => 4,
            );
        unless ($rc_cmd->execute) {
            die "Failed to obtain readcounts for file $file.\n";
        }
    } else {
        Genome::Sys->shellcmd( cmd => "touch $file.rcnt" );
    }
    return("$file.rcnt");
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
    unless (defined $model) {
        die "ERROR: Could not find a soamtic variation model with ID: " . $model->id . "\n";
    }
    unless (-e $output_dir) {
        die "ERROR: Output directory not found: $output_dir\n";
    }

    #grab the info from the model
    my $build = $model->last_succeeded_build;
    unless (defined($build)) {
        die "ERROR: Model ", $model->id, "has no succeeded builds\n";
    }

    my $tumor_model = $model->tumor_model;
    my $normal_model = $model->normal_model;

    my $ref_seq_build = $tumor_model->reference_sequence_build;
    my $ref_seq_fasta = $ref_seq_build->full_consensus_path('fa');
    my $annotation_build_name = $model->annotation_build->name;
    my $tiering_files = $model->annotation_build->data_directory . "/annotation_data/tiering_bed_files_v3/";
    my $sample_name = $self->sample_name;
    unless (defined($sample_name)) {
        $sample_name = $model->subject->name;
    }
    print STDERR "processing model with sample_name: " . $sample_name . "\n";
    my $tumor_bam = $tumor_model->last_succeeded_build->merged_alignment_result->bam_file;
    my $normal_bam = $normal_model->last_succeeded_build->merged_alignment_result->bam_file;
    my $build_dir = $build->data_directory;

    my $igv_reference_name = $self->igv_reference_name;
    if ($self->create_review_files && !defined($self->igv_reference_name)) {
        die "ERROR: igv-reference-name required if --create-review-files is specified";
    }


    # create subdirectories, get files in place

    # if multiple models with the same name, add a suffix
    if( -e "$output_dir/$sample_name" ){
        my $suffix = 1;
        my $newname = $sample_name . "-" . $suffix;
        while ( -e "$output_dir/$newname" ){
            $suffix++;
            $newname = $sample_name . "-" . $suffix;
        }
        $sample_name = $newname;
    }
    #make the directory structure
    Genome::Sys->create_directory("$output_dir/$sample_name");
    Genome::Sys->create_directory("$output_dir/$sample_name/snvs");
    Genome::Sys->create_directory("$output_dir/$sample_name/indels");
    unless ($self->create_review_files) {
        Genome::Sys->create_directory("$output_dir/review");
    }
    Genome::Sys->create_symlink($build_dir, "$output_dir/$sample_name/build_directory");


    # Check if the necessary files exist in this build
    my $snv_file = "$build_dir/effects/snvs.hq.novel.tier1.v2.bed";
    unless( -e $snv_file ){
        die "ERROR: SNV results for $sample_name not found at $snv_file\n";
    }
    my $indel_file = "$build_dir/effects/indels.hq.novel.tier1.v2.bed";
    unless( -e $indel_file ){
        die "ERROR: INDEL results for $sample_name not found at $indel_file\n";
    }
    my $sv_file;
    my $process_svs = $self->process_svs;
    if($process_svs){
        my @sv_files = glob("$build_dir/variants/sv/union-union-sv_breakdancer_*sv_squaredancer*/svs.merge.file.somatic");
        $sv_file = $sv_files[0];
        unless( -e $sv_file ){
            print STDERR "WARNING: SV results for $sample_name not found, skipping SVs\n";
#SHOULD THIS BE MADE INTO AN OBJECT PARAM $self->process_svs
            $process_svs = 0;
        }
    }


    #cat all the filtered snvs together (same for indels)
    Genome::Sys->shellcmd( cmd => "cat $build_dir/effects/snvs.hq.novel.tier*.v2.bed $build_dir/effects/snvs.hq.previously_detected.tier*.v2.bed | joinx sort >$output_dir/$sample_name/snvs/snvs.hq.bed");
    Genome::Sys->shellcmd( cmd => "cat $build_dir/effects/indels.hq.novel.tier*.v2.bed $build_dir/effects/indels.hq.previously_detected.tier*.v2.bed | joinx sort >$output_dir/$sample_name/indels/indels.hq.bed" );

#  `ln -s $snv_file $output_dir/$sample_name/snvs/` unless( -e "$output_dir/$sample_name/snvs/$snv_file");
#  `ln -s $indel_file $output_dir/$sample_name/indels/` unless( -e "$output_dir/$sample_name/indels/$indel_file");
  if($process_svs){
      Genome::Sys->create_directory("$output_dir/$sample_name/svs");
      Genome::Sys->create_symlink($sv_file, "$output_dir/$sample_name/svs/svs.hq") unless( -e "$output_dir/$sample_name/svs/$sv_file");

#       #annotate the svs
#       my $anno_cmd = Genome::Model::Tools::Annotate::Sv::Combine->create(
#           input-file => $sv_file,
#           output-file => "$output_dir/$sample_name/svs/svs.hq.annotated",
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
  $snv_file = "$output_dir/$sample_name/snvs/snvs.hq.bed";
  $indel_file = "$output_dir/$sample_name/indels/indels.hq.bed";



  my %dups;


  #munge through SNV file to remove duplicates and fix IUB codes
  #--------------------------------------------------------------
  $snv_file = cleanFile($snv_file);
  $indel_file = cleanFile($indel_file);



  #-------------------------------------------------
  #filter out the off-target regions, if target regions are available
  if($self->restrict_to_target_regions){
      print STDERR "Filtering out off-target regions...\n";

      my $featurelist;
      if($self->target_regions) {
         $featurelist = $self->target_regions;
      } elsif ($model->tumor_model->can('target_region_set_name') and defined($model->tumor_model->target_region_set_name)){
          my $featurelist_name = $model->tumor_model->target_region_set_name;
          $featurelist = Genome::FeatureList->get(name=>$featurelist_name)->file_path;
      }
      if(defined($featurelist) && (-s $featurelist)){
          #clean up feature list
          open(FEATFILE,">$output_dir/$sample_name/featurelist.tmp");
          my $inFh = IO::File->new( $featurelist ) || die "can't open file feature file\n";
          while( my $line = $inFh->getline )
          {
              chomp($line);
              next if $line =~ /^track/;
              my ( $chr, $start, $stop, @rest) = split( /\t/, $line );
              #remove chr if present
              $chr =~ s/^chr//g;
              print FEATFILE join("\t",( $chr, $start, $stop, @rest)) . "\n";
          }
          close($inFh);
          close(FEATFILE);
          my $new_snv_file = addName($snv_file,"ontarget");
          my $new_indel_file = addName($indel_file,"ontarget");

          Genome::Sys->shellcmd(cmd => "joinx sort $output_dir/$sample_name/featurelist.tmp >$output_dir/$sample_name/featurelist");
          Genome::Sys->shellcmd(cmd => "rm -f $output_dir/$sample_name/featurelist.tmp");
          Genome::Sys->shellcmd(cmd => "joinx intersect -a $snv_file -b $output_dir/$sample_name/featurelist >$new_snv_file");
          $snv_file = "$new_snv_file";
          Genome::Sys->shellcmd(cmd => "joinx intersect -a $indel_file -b $output_dir/$sample_name/featurelist >$new_indel_file");
          $indel_file = "$new_indel_file";
      }
      else {
          $self->warning_message("feature list not found or target regions not specified; No target region filtering being done even though --restrict-to-target-regions set.");
      }
  }

  #-------------------------------------------------
  #remove filter sites specified by the user
  if(defined($self->filter_sites)){
      print STDERR "Applying user-supplied filter...\n";
      my $filterSites = getFilterSites($self->filter_sites);
      $snv_file = removeFilterSites($snv_file,$filterSites);
      $indel_file = removeFilterSites($indel_file,$filterSites);
  }


    #-------------------------------------------------
    #remove filter regions specified by the user
    if(defined($self->filter_regions)){
        $snv_file = $self->_filter_regions($snv_file);
        $indel_file = $self->_filter_regions($indel_file);
    }


  #-------------------------------------------------------
  # remove regions called by less than the required number of callers
  unless($self->required_snv_callers == 1){
      print STDERR "Removing regions supported by less than " . $self->required_snv_callers . " regions...\n";
      $snv_file = removeUnsupportedSites($snv_file, $self->required_snv_callers, $build_dir);
  }




  ##------------------------------------------------------
  # do annotation
  $snv_file = doAnnotation($snv_file, $annotation_build_name);
  $indel_file = doAnnotation($indel_file, $annotation_build_name);


  #-------------------------------------------------------
  #add tiers
  if($self->add_tiers){
      print STDERR "Adding tiers...\n";
      #do annotation
      $snv_file = addTiering($snv_file, $tiering_files);
      $indel_file = addTiering($indel_file, $tiering_files);
  }



    #----------------------------------------------------
    # add dbsnp/gmaf
    if ($self->add_dbsnp_and_gmaf){
        $snv_file = $self->_add_dbsnp_and_gmaf_to_snv($build_dir, $snv_file);
        $indel_file = $self->_add_dbsnp_and_gmaf_to_indel($build_dir, $indel_file);
    }



  #-------------------------------------------------------
  #get readcounts
  if($self->get_readcounts){
      print STDERR "Getting readcounts...\n";
      if( -s "$snv_file" ){
          $snv_file = getReadcounts($snv_file, $ref_seq_fasta, $normal_bam, $tumor_bam);
      }
      if( -s "$indel_file" ){
          $indel_file = getReadcounts($indel_file, $ref_seq_fasta, $normal_bam, $tumor_bam);
      }
  }



    #------------------------------------------------------
    # combine the files into one master table
    $self->_create_master_files($snv_file, $indel_file);

    #------------------------------------------------------
    # now get the files together for review
    if($self->create_review_files){
        $self->_create_review_files();
    }

    #------------------------------------------------
    # tar up the files to be sent to collaborators
    my $archive_dir = "$output_dir/$sample_name/$sample_name";
    if($self->create_archive){
        $self->_create_archive($archive_dir, $build_dir, $sv_file);
    }

    return 1;
}

sub _filter_regions {
    my $self = shift;
    my $file = shift;

    #cat all the filter files together into one bed file
    #This will be executed twice because this subroutine is being called twice
    #which is not ideal
    my @filters = split(",",$self->filter_regions);

    my ($fh,$temp_file) = Genome::Sys->create_temp_file;

    foreach my $filterfile (@filters) {
        my $inFh2 = IO::File->new( $filterfile ) || die "can't open file8\n";
        while ( my $line = $inFh2->getline ) {
            print $fh $line;
        }
        close $inFh2;
    }
    close($fh);

    print STDERR "Removing user-specified filter for $file...\n";
    my $cmd = "joinx intersect --miss-a $file.filteredReg -a $file -b $temp_file >/dev/null";
    my $result = Genome::Sys->shellcmd(
        cmd => "$cmd",
    );
    unless ($result) {
        $self->error_message("Failed to execute joinx: Returned $result");
        die $self->error_message;
    }
    $file = "$file.filteredReg";
    return $file;
}

sub _add_dbsnp_and_gmaf_to_snv {
    my $self       = shift;
    my $build_dir  = shift;
    my $snv_file   = shift;

    print STDERR "==== adding dbsnp ids ====\n";
    print STDERR "$build_dir/variants/snvs.annotated.vcf.gz\n";
    if (-s "$build_dir/variants/snvs.annotated.vcf.gz") {
        my $db_cmd = Genome::Model::Tools::Annotate::AddRsid->create(
            anno_file => $snv_file,
            output_file => "$snv_file.rsid",
            vcf_file => "$build_dir/variants/snvs.annotated.vcf.gz",
        );
        unless ($db_cmd->execute) {
            die "Failed to add dbsnp anno to file $snv_file.\n";
        }
        $snv_file = "$snv_file.rsid";
        return $snv_file;
    }
    else {
        print STDERR "Warning: couldn't find annotated SNV file in build, skipping dbsnp anno\n";
    }
}

sub _add_dbsnp_and_gmaf_to_indel {
    my $self       = shift;
    my $build_dir  = shift;
    my $indel_file = shift;

    print STDERR "==== padding indel file with tabs to match ====\n";
    print STDERR "$build_dir/variants/snvs.annotated.vcf.gz\n";
    if (-s "$build_dir/variants/snvs.annotated.vcf.gz") {
        #pad indel file with tabs to match - if we ever start annotating with indels from dbsnp, replace this section
        open(OUTFILE,">$indel_file.rsid");
        my $inFh = IO::File->new( "$indel_file" ) || die "can't open file2d\n";
        while ( my $line = $inFh->getline ) {
            chomp($line);
            print OUTFILE $line . "\t\t\n"
        }
        close($inFh);
        close(OUTFILE);

        $indel_file = "$indel_file.rsid";
        return $indel_file;
    }
    else {
        print STDERR "Warning: couldn't find annotated SNV file in build, skipping dbsnp anno\n";
    }
}

sub _create_master_files {
    my $self       = shift;
    my $snv_file   = shift;
    my $indel_file = shift;

    my $output_dir = $self->output_dir;
    my $sample_name = $self->sample_name;
    unless (defined($sample_name)) {
        $sample_name = $self->somatic_variation_model->subject->name;
    }
    my $full_output_dir = "$output_dir/$sample_name";

    Genome::Sys->shellcmd(cmd => "head -n 1 $snv_file >$full_output_dir/snvs.indels.annotated");
    Genome::Sys->shellcmd(cmd => "tail -n +2 $indel_file >>$full_output_dir/snvs.indels.annotated.tmp");
    Genome::Sys->shellcmd(cmd => "tail -n +2 $snv_file >>$full_output_dir/snvs.indels.annotated.tmp");
    Genome::Sys->shellcmd(cmd => "joinx sort -i $full_output_dir/snvs.indels.annotated.tmp >>$full_output_dir/snvs.indels.annotated");
    Genome::Sys->shellcmd(cmd => "rm -f $full_output_dir/snvs.indels.annotated.tmp");

    # convert master table to excel
    my $workbook  = Spreadsheet::WriteExcel->new("$output_dir/$sample_name/snvs.indels.annotated.xls");
    my $worksheet = $workbook->add_worksheet();

    my $row=0;
    my $inFh = IO::File->new( "$output_dir/$sample_name/snvs.indels.annotated" ) || die "can't open file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my @F = split("\t",$line);
        for(my $i=0;$i<@F;$i++){
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

    my $output_dir = $self->output_dir;
    my $sample_name = $self->sample_name;
    unless (defined($sample_name)) {
        $sample_name = $self->somatic_variation_model->subject->name;
    }
    my $full_output_dir = "$output_dir/$sample_name";

    print "Generating Review files...\n";
    my @tiers = split(",",$self->tiers_to_review);
    my $tierstring = join("",@tiers);
    for my $i (@tiers){
        Genome::Sys->shellcmd(cmd => "grep -w tier$i $full_output_dir/snvs.indels.annotated >>$full_output_dir/snvs.indels.annotated.tier$tierstring.tmp");
    }
    Genome::Sys->shellcmd(cmd => "joinx sort -i $full_output_dir/snvs.indels.annotated.tier$tierstring.tmp >$full_output_dir/snvs.indels.annotated.tier$tierstring");
    Genome::Sys->shellcmd(cmd => "rm -f $output_dir/$sample_name/snvs.indels.annotated.tier$tierstring.tmp");
    annoFileToSlashedBedFile("$full_output_dir/snvs.indels.annotated.tier$tierstring","$output_dir/review/$sample_name.bed");

    my $tumor_bam = $self->somatic_variation_model->tumor_model->last_succeeded_build->merged_alignment_result->bam_file;
    my $normal_bam = $self->somatic_variation_model->normal_model->last_succeeded_build->merged_alignment_result->bam_file;

    my $bam_files = join(",",($normal_bam,$tumor_bam));
    my $labels = join(",",("normal $sample_name","tumor $sample_name"));

    my $igv_reference_name = $self->igv_reference_name;
    unless (defined($igv_reference_name))  {
          print STDERR "WARNING: No IGV reference name supplied - defaulting to build 37\n";
    }

    #create the xml file for review
    my $dumpXML = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
        bams => "$bam_files",
        labels => "$labels",
        output_file => "$output_dir/review/$sample_name.xml",
        genome_name => $sample_name,
        review_bed_file => "$output_dir/review/$sample_name.bed",
        reference_name => $igv_reference_name,
    );
    unless ($dumpXML->execute) {
        die "Failed to create IGV xml file\n";
    }

    print STDERR "\n--------------------------------------------------------------------------------\n";
    print STDERR "Sites to review are here:\n";
    print STDERR "$output_dir/review/$sample_name.bed\n";
    print STDERR "IGV XML file is here:";
    print STDERR "$output_dir/review/$sample_name.xml\n\n";

    return 1;
}

sub _create_archive {
    my $self        = shift;
    my $archive_dir = shift;
    my $build_dir   = shift;
    my $sv_file     = shift;

    my $output_dir = $self->output_dir;

    my $sample_name = $self->sample_name;
    unless (defined($sample_name)) {
        $sample_name = $self->somatic_variation_model->subject->name;
    }

    my $full_output_dir = "$output_dir/$sample_name";

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
            print STDERR "WARNING: no indel VCF file available. If this is an older model, a rebuild may fix this\n";
        }

        if (-e "$build_dir/variants/snvs.annotated.vcf.gz") {
            Genome::Sys->create_symlink("$build_dir/variants/snvs.annotated.vcf.gz", "$archive_dir/snvs.vcf.gz");
        }
        elsif (-e "$build_dir/variants/snvs.vcf.gz") {
            Genome::Sys->create_symlink("$build_dir/variants/snvs.vcf.gz", "$archive_dir/snvs.vcf.gz");
        }
        else {
            print STDERR "WARNING: no snv VCF file available. If this is an older model, a rebuild may fix this\n";
        }
    }

    #symlink annotated snvs and indels
    Genome::Sys->create_symlink("$full_output_dir/snvs.indels.annotated", "$archive_dir/snvsAndIndels.annotated");
    #symlink annoted snvs and indels excel file
    Genome::Sys->create_symlink("$full_output_dir/snvs.indels.annotated.xls", "$archive_dir/snvsAndIndels.annotated.xls");
    #symlink sv calls
    if($self->process_svs){
        Genome::Sys->create_symlink("$sv_file", "$archive_dir/svs");
    }
    #synlink cnv calls - todo

    #tar it up
#WRITE Genome::Sys TEST FOR THE NEW options PARAMETER
    Genome::Sys->tar(tar_path => "$archive_dir.tar.gz", input_directory => "$archive_dir", options => "-chzvf");

    return 1;
}

1;
