package Genome::Model::ClinSeq::Command::GetVariantSources;

#Written by Obi Griffith

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::GetVariantSources {
    is        => 'Command::V2',
    has_input => [
        builds => {
            is                  => 'Genome::Model::Build::SomaticInterface',
            is_many             => 1,
            shell_args_position => 1,
            require_user_verify => 0,
            doc                 => 'somatic build(s) to get variant sources from',
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
    doc => 'summarize the sources of variants (i.e., which variant callers) for a somatic build',
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
Summarize source of variants (i.e., snv/indel caller) for one or more somatic builds

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

    for my $somatic_build ($self->builds) {
        my $build_outdir = $self->build_outdir($somatic_build);

        my $indel_results_file = $self->concatenate_indel_files($somatic_build);
        my $snv_results_file   = $self->concatenate_snv_files($somatic_build);

        #Sort the BED files using joinx
        my $indel_results_file_sorted = $indel_results_file . '.sort';
        $self->joinxSortFile($indel_results_file, $indel_results_file_sorted);
        unlink $indel_results_file;
        Genome::Sys->move_file($indel_results_file_sorted, $indel_results_file);

        my $snv_results_file_sorted = $snv_results_file . '.sort';
        $self->joinxSortFile($snv_results_file, $snv_results_file_sorted);
        unlink $snv_results_file;
        Genome::Sys->move_file($snv_results_file_sorted, $snv_results_file);

        #Locate the individual indel/snv files for each caller to use in joinx intersect
        #This should be replaced by a method which somehow determines the appropriate files automatically
        #Depending on whether the somatic build is for exome or wgs data, the paths will differ - this should also be determined automatically
        for my $caller (qw(strelka gatk pindel varscan)) {
            $self->get_variant_caller_results("indel", $caller, $somatic_build);
        }

        for my $caller (qw(strelka sniper varscan samtools mutect)) {
            $self->get_variant_caller_results("snv", $caller, $somatic_build);
        }

        my %indel_caller_stats = %{$self->createIndelOutfile($somatic_build)};
        my %snv_caller_stats   = %{$self->createSnvOutfile($somatic_build)};
        $self->writeStats($build_outdir, \%snv_caller_stats, \%indel_caller_stats);
        #Cleanup temp files
        unlink $indel_results_file;
        unlink $snv_results_file;
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

sub build_outdir {
    my $self          = shift;
    my $somatic_build = shift;

    #If there is more than one somatic build supplied... create sub-directories for each
    my $build_outdir;
    if (scalar($self->builds) > 1) {
        $build_outdir = File::Spec->join($self->outdir, $somatic_build->id);
        unless (-e $build_outdir) {
            Genome::Sys->create_directory($build_outdir);
        }
    }
    else {
        $build_outdir = $self->outdir;
    }

    return $build_outdir;
}

sub snv_results_file {
    my $self  = shift;
    my $build = shift;

    return File::Spec->join($self->build_outdir($build), 'snvs.hq.novel.tier1-3.v2.bed');
}

sub indel_results_file {
    my $self  = shift;
    my $build = shift;

    return File::Spec->join($self->build_outdir($build), 'indels.hq.novel.tier1-3.v2.bed');
}

sub create_snv_hash {
    my $self          = shift;
    my $somatic_build = shift;

    my %snvs = %{$self->create_snv_indel_hash($self->get_snv_files_names($somatic_build))};
    $self->debug_message("Stored %s indels", scalar(keys %snvs));

    return %snvs;
}

sub create_indel_hash {
    my $self         = shift;
    my $somatic_build = shift;

    my %indels = %{$self->create_snv_indel_hash($self->get_indel_files_names($somatic_build))};
    $self->debug_message("Stored %s indels", scalar(keys %indels));

    return %indels;
}

sub create_snv_indel_hash {
    my $self  = shift;
    my $files = shift;

    my $variants;
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
        }
        close(VARIANTS);
    }

    return $variants;
}

sub concatenate_snv_files {
    my $self        = shift;
    my $somatic_build = shift;

    my %files = %{$self->get_snv_files_names($somatic_build)};
    my $results_file = $self->snv_results_file($somatic_build);
    my @file_list;
    for my $c (sort {$a <=> $b} keys %files) {
        push @file_list, $files{$c}{file};
    }
    Genome::Sys->concatenate_files(\@file_list, $results_file);

    return $results_file;
}

sub concatenate_indel_files {
    my $self        = shift;
    my $somatic_build = shift;

    my %files = %{$self->get_indel_files_names($somatic_build)};
    my $results_file = $self->indel_results_file($somatic_build);
    my @file_list;
    for my $c (sort {$a <=> $b} keys %files) {
        push @file_list, $files{$c}{file};
    }
    Genome::Sys->concatenate_files(\@file_list, $results_file);\

    return $results_file;
}

sub createIndelOutfile {
    my $self               = shift;
    my $somatic_build      = shift;

    my %indels        = $self->create_indel_hash($somatic_build);
    my $indel_outfile = File::Spec->join($self->build_outdir($somatic_build), 'indel_sources.tsv');

    my $indel_caller_stats;
    foreach my $caller ("strelka", "gatk", "pindel", "varscan") {
        $indel_caller_stats->{$caller} = 0;
    }
    open(INDEL_OUT, ">$indel_outfile") || die "\n\nCould not open $indel_outfile\n\n";
    print INDEL_OUT "coord\tchr\tstart\tend\tvariant\tscore1\tscore2\tcallers\tstrelka\tgatk\tpindel\tvarscan\ttier\n";
    foreach my $indel (sort {$indels{$a}->{coord_string} cmp $indels{$b}->{coord_string}} keys %indels) {
        my @callers = sort keys %{$indel_caller{$indels{$indel}{variant_string}}};
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
        print INDEL_OUT "$indels{$indel}{coord_string}\t$indels{$indel}{line}\t", join(",", @callers),
            "\t$strelka\t$gatk\t$pindel\t$varscan\t$indels{$indel}{tier}\n";
    }
    close(INDEL_OUT);

    #Set output files as output to this step
    die $self->error_message("Trying to set a file as output but the file does not exist: $indel_outfile")
        unless (-e $indel_outfile);
    $self->indel_variant_sources_file($indel_outfile);

    return $indel_caller_stats;
}

sub createSnvOutfile {
    my $self             = shift;
    my $somatic_build    = shift;

    my %snvs          = $self->create_snv_hash($somatic_build);
    my $snv_outfile   = File::Spec->join($self->build_outdir($somatic_build), 'snv_sources.tsv');

    my $snv_caller_stats;
    foreach my $caller ("strelka", "sniper", "varscan", "samtools", "mutect") {
        $snv_caller_stats->{$caller} = 0;
    }
    open(SNV_OUT, ">$snv_outfile") || die "\n\nCould not open $snv_outfile\n\n";
    print SNV_OUT
        "coord\tchr\tstart\tend\tvariant\tscore1\tscore2\tcallers\tstrelka\tsniper\tvarscan\tsamtools\tmutect\ttier\n";
    foreach my $snv (sort {$snvs{$a}->{coord_string} cmp $snvs{$b}->{coord_string}} keys %snvs) {
        my @callers  = sort keys %{$snv_caller{$snvs{$snv}{variant_string}}};
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
        print SNV_OUT "$snvs{$snv}{coord_string}\t$snvs{$snv}{line}\t", join(",", @callers),
            "\t$strelka\t$sniper\t" . "$varscan\t$samtools\t$mutect\t$snvs{$snv}{tier}\n";
    }
    close(SNV_OUT);

    #Set output files as output to this step
    die $self->error_message("Trying to set a file as output but the file does not exist: $snv_outfile")
        unless (-e $snv_outfile);
    $self->snv_variant_sources_file($snv_outfile);

    return $snv_caller_stats;
}

sub get_indel_files_names {
    my $self        = shift;
    my $build       = shift;

    my $indel_files;
    for my $tier (1..4) {
        my $indel_file = $build->indels_effects_file("tier$tier");
        $indel_files->{$tier}{file} = $indel_file;
        $indel_files->{$tier}{tier} = "tier$tier";
    }

    return $indel_files;
}

sub get_snv_files_names {
    my $self      = shift;
    my $build     = shift;

    my $snv_files;
    for my $tier (1..4) {
        my $snv_file = $build->snvs_effects_file("tier$tier");
        $snv_files->{$tier}{file} = $snv_file;
        $snv_files->{$tier}{tier} = "tier$tier";
    }

    return $snv_files;
}

sub get_variant_caller_results {
    my $self                 = shift;
    my $variant_type         = shift;  #"snv" or "indel"
    my $caller               = shift;
    my $somatic_build        = shift;

    my $build_outdir         = $self->build_outdir($somatic_build);
    my $variant_caller_file;

    #Create a list of possible indels file paths
    my @variant_caller_paths =
        glob(File::Spec->join($somatic_build->data_directory, 'variants', $variant_type, "$caller-*", "${variant_type}s.hq.bed"));
    if (@variant_caller_paths) {
        $variant_caller_file = $self->checkResultFile(
            '-paths'  => \@variant_caller_paths,
            '-caller' => $caller
        );

        #Sort the caller result BED files using joinx and store in a temporary file and use that to run joinx intersect
        my $variant_caller_file_s = File::Spec->join($build_outdir, "${variant_type}_${caller}.sorted.bed");
        $self->joinxSortFile($variant_caller_file, $variant_caller_file_s);

        #Use 'joinx intersect' to determine which indels in the merged/union file are found in each individual caller's results file
        #gmt joinx intersect a.bed b.bed [--output-file=n.bed] --exact-pos --exact-allele
        my $params_string          = "--exact-pos --exact-allele";
        my $variant_caller_outfile = File::Spec->join($build_outdir, "${variant_type}_${caller}.bed");

        my $results_file_accessor = "${variant_type}_results_file";
        my $variant_results_file  = $self->$results_file_accessor($somatic_build);
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
