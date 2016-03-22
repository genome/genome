package Genome::Model::ClinSeq::Command::GetVariantSources;

#Written by Obi Griffith

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::GetVariantSources {
    is        => 'Command::V2',
    has_input => [
        builds => {
            is                  => 'Genome::Model::Build::SomaticVariation',
            is_many             => 1,
            shell_args_position => 1,
            require_user_verify => 0,
            doc                 => 'somatic variation build(s) to get variant sources from',
        },
        outdir => {
            is  => 'FilesystemPath',
            doc => 'Directory where output files will be written',
        },
    ],
    has_output => [
        indel_variant_sources_file => {
            is          => 'FilesystemPath',
            is_optional => 1,
        },
        snv_variant_sources_file => {
            is          => 'FilesystemPath',
            is_optional => 1,
        },
    ],
    doc => 'summarize the sources of variants (i.e., which variant callers) for a somatic variation build',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq get-variant-sources --outdir=/tmp/  128884819

genome model clin-seq get-variant-sources --outdir=/tmp/  id=128884819

genome model clin-seq get-variant-sources --outdir=/tmp/  model.id=2888329352

genome model clin-seq get-variant-sources --outdir=/tmp/  "model.name='H_JG-300000-1206887.somatic_variation-1'"

genome model clin-seq get-variant-sources --outdir=/tmp/  'id in [128884819,128884852]'

EOS
}

sub help_detail {
    return <<EOS
Summarize source of variants (i.e., snv/indel caller) for one or more somatic variation builds

(put more content here)
EOS
}

sub __errors__ {
    my $self   = shift;
    my @errors = $self->SUPER::__errors__(@_);

    unless (-e $self->outdir && -d $self->outdir) {
        push @errors,
            UR::Object::Tag->create(
            type       => 'error',
            properties => ['outdir'],
            desc       => "Outdir: " . $self->outdir . " not found or not a directory",
            );
    }
    return @errors;
}

#Global variables - #figure out a way to do this with a method
my %indel_caller;
my %snv_caller;

sub execute {
    my $self   = shift;
    my @builds = $self->builds;
    my $outdir = $self->outdir;

    unless ($outdir =~ /\/$/) {
        $outdir .= "/";
    }

    my $somatic_build_count = scalar(@builds);
    for my $somatic_build (@builds) {
        #If there is more than one somatic variation build supplied... create sub-directories for each
        my $build_outdir;
        if ($somatic_build_count > 1) {
            $build_outdir = $outdir . $somatic_build->id . "/";
            mkdir($build_outdir);
        }
        else {
            $build_outdir = $outdir;
        }
        my $somatic_build_dir = $somatic_build->data_directory;

        #Set files for output
        my $indel_outfile = $build_outdir . "indel_sources.tsv";
        my $snv_outfile   = $build_outdir . "snv_sources.tsv";

        #Set a list of variant files to consider
        my %indel_files;
        $self->get_indel_files_names(\%indel_files, $somatic_build_dir);

        my %snv_files;
        $self->get_snv_files_names(\%snv_files, $somatic_build_dir);

        #Locate the final indel/snv results files and load into memory
        #For indels, use ~/effects/indels.hq.novel.tier1.v2.bed ?  (Or the annotated file?)
        #For SNVs, use ~/effects/snvs.hq.novel.tier1.v2.bed
        my (%indels, %snvs);
        my $indel_results_file = $build_outdir . "indels.hq.novel.tier1-3.v2.bed";
        my $snv_results_file   = $build_outdir . "snvs.hq.novel.tier1-3.v2.bed";
        $self->create_snv_indel_hash(\%indel_files, \%indels, $indel_results_file);
        $self->create_snv_indel_hash(\%snv_files,   \%snvs,   $snv_results_file);

        #Sort the BED files using joinx
        my $indel_results_file_sorted = $indel_results_file . ".sort";
        $self->joinxSortFile($indel_results_file, $indel_results_file_sorted);
        unlink $indel_results_file;
        Genome::Sys->move_file($indel_results_file_sorted, $indel_results_file);

        my $snv_results_file_sorted = $snv_results_file . ".sort";
        $self->joinxSortFile($snv_results_file, $snv_results_file_sorted);
        unlink $snv_results_file;
        Genome::Sys->move_file($snv_results_file_sorted, $snv_results_file);

        my $indel_count = keys %indels;
        $self->debug_message("Stored $indel_count indels");
        my $snv_count = keys %snvs;
        $self->debug_message("Stored $snv_count indels");

        #Locate the individual indel/snv files for each caller to use in joinx intersect
        #This should be replaced by a method which somehow determines the appropriate files automatically
        #Depending on whether the somatic variation build is for exome or wgs data, the paths will differ - this should also be determined automatically
        $self->get_variant_caller_results("indel", "strelka", $indel_results_file, $somatic_build_dir, $build_outdir);
        $self->get_variant_caller_results("indel", "gatk",    $indel_results_file, $somatic_build_dir, $build_outdir);
        $self->get_variant_caller_results("indel", "pindel",  $indel_results_file, $somatic_build_dir, $build_outdir);
        $self->get_variant_caller_results("indel", "varscan", $indel_results_file, $somatic_build_dir, $build_outdir);

        $self->get_variant_caller_results("snv", "strelka",  $snv_results_file, $somatic_build_dir, $build_outdir);
        $self->get_variant_caller_results("snv", "sniper",   $snv_results_file, $somatic_build_dir, $build_outdir);
        $self->get_variant_caller_results("snv", "varscan",  $snv_results_file, $somatic_build_dir, $build_outdir);
        $self->get_variant_caller_results("snv", "samtools", $snv_results_file, $somatic_build_dir, $build_outdir);
        $self->get_variant_caller_results("snv", "mutect",   $snv_results_file, $somatic_build_dir, $build_outdir);

        my (%indel_caller_stats, %snv_caller_stats);
        $self->createIndelOutfile(\%indels, $indel_outfile, \%indel_caller_stats);
        $self->createSnvOutfile(\%snvs, $snv_outfile, \%snv_caller_stats);
        $self->writeStats($build_outdir, \%snv_caller_stats, \%indel_caller_stats);
        #Cleanup temp files
        unlink $indel_results_file;
        unlink $snv_results_file;

        #Set output files as output to this step
        die $self->error_message("Trying to set a file as output but the file does not exist: $indel_outfile")
            unless (-e $indel_outfile);
        $self->indel_variant_sources_file($indel_outfile);
        die $self->error_message("Trying to set a file as output but the file does not exist: $snv_outfile")
            unless (-e $snv_outfile);
        $self->snv_variant_sources_file($snv_outfile);
    }

    $self->debug_message("\n\n");

    return 1;
}

sub noteCaller {
    my $self           = shift;
    my $intersect_file = shift;
    my $caller         = shift;
    my $variant_type   = shift;
    #Go through output of 'joinx intersect' and identify which lines in the main/original file were found intersecting with the other file
    open(INTERSECT, "$intersect_file") || die "\n\ncan't open $intersect_file\n";

    while (<INTERSECT>) {
        chomp;
        my @data = split("\t", $_);
        my $variant_string = "$data[0]" . ":" . "$data[1]" . "-" . "$data[2]" . " ($data[3])";
        if ($variant_type eq "indel") {
            $indel_caller{$variant_string}{$caller}++;
        }
        if ($variant_type eq "snv") {
            $snv_caller{$variant_string}{$caller}++;
        }
    }
    #when finished, delete intermediate result file from joinx - otherwise this can sometimes interfere with future runs of the tool
    unlink $intersect_file;
}

sub checkResultFile {
    my $self        = shift;
    my %args        = @_;
    my @paths       = @{$args{'-paths'}};
    my $caller      = $args{'-caller'};
    my $result_file = '';
    foreach my $path (@paths) {
        if (-e $path) {
            $result_file = $path;
        }
    }
    unless (-e $result_file) {
        my $path_list = join("\n", @paths);
        $self->error_message("$caller result not found in the following list of paths\n\n$path_list\n");
        return undef;
    }
    return ($result_file);
}

sub create_snv_indel_hash {
    my $self         = shift;
    my $files        = shift;
    my $variants     = shift;
    my $results_file = shift;
    open(OUT, ">$results_file") || die $self->error_message("Could not open output file: $results_file");
    my $l = 0;
    foreach my $c (sort {$a <=> $b} keys %$files) {
        my $file = $files->{$c}{file};
        my $tier = $files->{$c}{tier};
        open(VARIANTS, $file) or die "can't open $file\n";
        while (<VARIANTS>) {
            $l++;
            chomp;
            my $line           = $_;
            my @data           = split("\t", $_);
            my $variant_string = "$data[0]" . ":" . "$data[1]" . "-" . "$data[2]" . " ($data[3])";
            my $coord_string   = "$data[0]" . ":" . "$data[1]" . "-" . "$data[2]";
            $variants->{$l}{line}           = $line;
            $variants->{$l}{variant_string} = $variant_string;
            $variants->{$l}{coord_string}   = $coord_string;
            $variants->{$l}{tier}           = $tier;
            print OUT "$line\n";
        }
        close(VARIANTS);
    }
    close(OUT);
}

sub createIndelOutfile {
    my $self               = shift;
    my $indels             = shift;
    my $indel_outfile      = shift;
    my $indel_caller_stats = shift;
    foreach my $caller ("strelka", "gatk", "pindel", "varscan") {
        $indel_caller_stats->{$caller} = 0;
    }
    open(INDEL_OUT, ">$indel_outfile") || die "\n\nCould not open $indel_outfile\n\n";
    print INDEL_OUT "coord\tchr\tstart\tend\tvariant\tscore1\tscore2\tcallers\tstrelka\tgatk\tpindel\tvarscan\ttier\n";
    foreach my $indel (sort {$indels->{$a}->{coord_string} cmp $indels->{$b}->{coord_string}} keys %$indels) {
        my @callers = sort keys %{$indel_caller{$indels->{$indel}{variant_string}}};
        my $strelka = 0;
        my $gatk    = 0;
        my $pindel  = 0;
        my $varscan = 0;
        foreach my $caller (@callers) {
            if ($caller eq 'strelka') {$strelka = 1;}
            if ($caller eq 'gatk')    {$gatk    = 1;}
            if ($caller eq 'pindel')  {$pindel  = 1;}
            if ($caller eq 'varscan') {$varscan = 1;}
            $indel_caller_stats->{$caller} += 1;
        }
        print INDEL_OUT "$indels->{$indel}{coord_string}\t$indels->{$indel}{line}\t", join(",", @callers),
            "\t$strelka\t$gatk\t$pindel\t$varscan\t$indels->{$indel}{tier}\n";
    }
    close(INDEL_OUT);
}

sub createSnvOutfile {
    my $self             = shift;
    my $snvs             = shift;
    my $snv_outfile      = shift;
    my $snv_caller_stats = shift;
    foreach my $caller ("strelka", "sniper", "varscan", "samtools", "mutect") {
        $snv_caller_stats->{$caller} = 0;
    }
    open(SNV_OUT, ">$snv_outfile") || die "\n\nCould not open $snv_outfile\n\n";
    print SNV_OUT
        "coord\tchr\tstart\tend\tvariant\tscore1\tscore2\tcallers\tstrelka\tsniper\tvarscan\tsamtools\tmutect\ttier\n";
    foreach my $snv (sort {$snvs->{$a}->{coord_string} cmp $snvs->{$b}->{coord_string}} keys %$snvs) {
        my @callers  = sort keys %{$snv_caller{$snvs->{$snv}{variant_string}}};
        my $strelka  = 0;
        my $sniper   = 0;
        my $varscan  = 0;
        my $samtools = 0;
        my $mutect   = 0;
        foreach my $caller (@callers) {
            if ($caller eq 'strelka')  {$strelka  = 1;}
            if ($caller eq 'sniper')   {$sniper   = 1;}
            if ($caller eq 'varscan')  {$varscan  = 1;}
            if ($caller eq 'samtools') {$samtools = 1;}
            if ($caller eq 'mutect')   {$mutect   = 1;}
            $snv_caller_stats->{$caller} += 1;
        }
        print SNV_OUT "$snvs->{$snv}{coord_string}\t$snvs->{$snv}{line}\t", join(",", @callers),
            "\t$strelka\t$sniper\t" . "$varscan\t$samtools\t$mutect\t$snvs->{$snv}{tier}\n";
    }
    close(SNV_OUT);
}

sub get_indel_files_names {
    my $self        = shift;
    my $indel_files = shift;
    my $build_dir   = shift;
    $indel_files->{1}{file} = $build_dir . "/effects/indels.hq.novel.tier1.v2.bed";
    $indel_files->{1}{tier} = "tier1";
    $indel_files->{2}{file} = $build_dir . "/effects/indels.hq.novel.tier2.v2.bed";
    $indel_files->{2}{tier} = "tier2";
    $indel_files->{3}{file} = $build_dir . "/effects/indels.hq.novel.tier3.v2.bed";
    $indel_files->{3}{tier} = "tier3";
    $indel_files->{4}{file} = $build_dir . "/effects/indels.hq.novel.tier4.v2.bed";
    $indel_files->{4}{tier} = "tier4";
}

sub get_snv_files_names {
    my $self      = shift;
    my $snv_files = shift;
    my $build_dir = shift;
    $snv_files->{1}{file} = $build_dir . "/effects/snvs.hq.novel.tier1.v2.bed";
    $snv_files->{1}{tier} = "tier1";
    $snv_files->{2}{file} = $build_dir . "/effects/snvs.hq.novel.tier2.v2.bed";
    $snv_files->{2}{tier} = "tier2";
    $snv_files->{3}{file} = $build_dir . "/effects/snvs.hq.novel.tier3.v2.bed";
    $snv_files->{3}{tier} = "tier3";
    $snv_files->{4}{file} = $build_dir . "/effects/snvs.hq.novel.tier4.v2.bed";
    $snv_files->{4}{tier} = "tier4";
}

sub get_variant_caller_results() {
    my $self                 = shift;
    my $variant_type         = shift;  #"snv" or "indel"
    my $caller               = shift;
    my $variant_results_file = shift;
    my $build_dir            = shift;
    my $build_outdir         = shift;
    my $variant_caller_file;

    #Create a list of possible indels file paths
    my @variant_caller_paths =
        glob("${build_dir}/variants/" . $variant_type . "/" . $caller . "-*/" . $variant_type . "s.hq.bed");
    if (@variant_caller_paths) {
        $variant_caller_file = $self->checkResultFile(
            '-paths'  => \@variant_caller_paths,
            '-caller' => $caller
        );

        #Sort the caller result BED files using joinx and store in a temporary file and use that to run joinx intersect
        my $variant_caller_file_s = $build_outdir . $variant_type . "_" . $caller . ".sorted.bed";
        $self->joinxSortFile($variant_caller_file, $variant_caller_file_s);

        #Use 'joinx intersect' to determine which indels in the merged/union file are found in each individual caller's results file
        #gmt joinx intersect a.bed b.bed [--output-file=n.bed] --exact-pos --exact-allele
        my $params_string          = "--exact-pos --exact-allele";
        my $variant_caller_outfile = $build_outdir . $variant_type . "_" . $caller . ".bed";

        $self->determineCaller($variant_type, $caller, $variant_results_file, $variant_caller_file_s,
            $variant_caller_outfile);
        unlink $variant_caller_file_s;
    }
}

sub determineCaller {
    my ($self, $variant_type, $caller_name, $results_file, $caller_file, $outfile) = @_;

    if (defined $caller_file) {
        $self->debug_message("Looking for overlapping $variant_type results between:\n$results_file\n$caller_file\n\n");

        my $cmd = Genome::Model::Tools::Joinx::Intersect->create(
            exact_pos    => 1,
            exact_allele => 1,
            output_file  => $outfile,
            input_file_a => $results_file,
            input_file_b => $caller_file
        );
        $cmd->execute();

        #Go through original indels and note all files from different callers where that indel was called
        $self->noteCaller($outfile, $caller_name, $variant_type);
    }
}

sub joinxSortFile {
    my ($self, $ip_file, $sorted_op_file) = @_;
    my $joinx_sort_cmd = Genome::Model::Tools::Joinx::Sort->create(
        output_file => $sorted_op_file,
        input_files => [$ip_file]
    );
    $joinx_sort_cmd->execute();
}

sub writeStats {
    my $self               = shift;
    my $outdir             = shift;
    my $snv_caller_stats   = shift;
    my $indel_caller_stats = shift;
    my $stats_file         = $outdir . "/Stats.tsv";
    my $data_source        = "WGS/Exome";
    if ($stats_file =~ /variant_source_callers\/exome/) {
        $data_source = "Exome";
    }
    elsif ($stats_file =~ /variant_source_callers\/wgs/) {
        $data_source = "WGS";
    }
    open my $STATS, ">$stats_file";
    print $STATS "Question\tAnswer\tData_Type\tAnalysis_Type\tStatistic_Type\tExtra_Description\n";
    $self->write_snv_stats($snv_caller_stats, $data_source, $STATS);
    $self->write_indel_stats($indel_caller_stats, $data_source, $STATS);
    close $STATS;
}

sub write_snv_stats {
    my $self             = shift;
    my $snv_caller_stats = shift;
    my $data_source      = shift;
    my $STATS            = shift;
    print $STATS "Number of Strelka SNV calls\t"
        . $snv_caller_stats->{"strelka"} . "\t"
        . $data_source
        . "\tClinseq Build Summary\tCount\tNumber of SNVs called by Strelka\n";
    print $STATS "Number of Sniper SNV calls\t"
        . $snv_caller_stats->{"sniper"} . "\t"
        . $data_source
        . "\tClinseq Build Summary\tCount\tNumber of SNVs called by Sniper\n";
    print $STATS "Number of VarScan SNV calls\t"
        . $snv_caller_stats->{"varscan"} . "\t"
        . $data_source
        . "\tClinseq Build Summary\tCount\tNumber of SNVs called by VarScan\n";
    print $STATS "Number of SamTools SNV calls\t"
        . $snv_caller_stats->{"samtools"} . "\t"
        . $data_source
        . "\tClinseq Build Summary\tCount\tNumber of SNVs called by SamTools\n";
    print $STATS "Number of Mutect SNV calls\t"
        . $snv_caller_stats->{"mutect"} . "\t"
        . $data_source
        . "\tClinseq Build Summary\tCount\tNumber of SNVs called by Mutect\n";
}

sub write_indel_stats {
    my $self               = shift;
    my $indel_caller_stats = shift;
    my $data_source        = shift;
    my $STATS              = shift;
    print $STATS "Number of Strelka Indel calls\t"
        . $indel_caller_stats->{"strelka"} . "\t"
        . $data_source
        . "\tClinseq Build Summary\tCount\tNumber of Indels called by Strelka\n";
    print $STATS "Number of GATK Indel calls\t"
        . $indel_caller_stats->{"gatk"} . "\t"
        . $data_source
        . "\tClinseq Build Summary\tCount\tNumber of Indels called by GATK\n";
    print $STATS "Number of Pindel Indel calls\t"
        . $indel_caller_stats->{"pindel"} . "\t"
        . $data_source
        . "\tClinseq Build Summary\tCount\tNumber of Indels called by Pindel\n";
    print $STATS "Number of VarScan Indel calls\t"
        . $indel_caller_stats->{"varscan"} . "\t"
        . $data_source
        . "\tClinseq Build Summary\tCount\tNumber of Indels called by Varscan\n";
}

1;
