package Genome::Model::Tools::Vcf::SomaticVariationToVcf;

use strict;
use warnings;
use Genome;
use File::stat;
use IO::File;
use File::Basename;
use Getopt::Long;
use FileHandle;
use Data::Dumper;

my %detector_to_module =
( samtools => 'vcf-maker-sniper-somatic',
    sniper   => 'vcf-maker-sniper-somatic',
    varscan  => 'vcf-maker-varscan-somatic',
    );
 $detector_to_module{'varscan-somatic'}='vcf-maker-varscan-somatic';
my %detector_vcf;

class Genome::Model::Tools::Vcf::SomaticVariationToVcf {
    is => 'Command',
    has => [
    output_dir => {
        is => 'Text',
        is_output => 1,
        is_optional => 0,
        doc => "Output merged VCF",
    },
    build_id => {
        is => 'Text',
        is_input=>1,
    },
    annotate_dbsnp => { default=>0 },
    genome_build => { default=>36 },
    ],
};


sub help_brief {                            # keep this to just a few words <---
    "Merge multiple VCFs - keep the quality scores from files in desc order"
}


sub help_synopsis {
    <<'HELP';
Merge multiple VCFs - keep the FORMAT lines from files in desc order.
HELP
}

sub help_detail {                  # this is what the user will see with the longer version of help. <---
    <<'HELP';
Merge multiple VCFs. For identical calls made by different algorithms, merge them, keeping the FORMAT/scores from the file that is listed first in the vcf_files string.
HELP
}

###############

sub execute {                               # replace with real execution logic
    my $self = shift;
    my $build = Genome::Model::Build->get($self->build_id);
    unless($build) {
        $self->error_message("Could not get object for build_id: " . $self->build_id);
        return;
    }

    my $dv2=Genome::Model::Tools::DetectVariants2::Dispatcher->create(snv_detection_strategy=>$build->snv_detection_strategy);
    my ($tree, $plan) = $dv2->plan;
    if($self->check_for_conversion_ability($plan)) {
        $self->status_message("All SNV detectors have a module to convert to vcf. Beginning conversion");
    }
    else {
         return 0;
    }
    my $build_dir = $build->data_directory;
    my $sample = $build->subject_name;

    ### for each detector, create the vcf then grab the results of each filter and apply them to that file
    for my $detector (keys %{$plan}) {
        my $module_name = $detector_to_module{$detector};
        my $output_file = $self->output_dir .  "/$sample.$detector.vcf";
        for my $version (keys %{$plan->{$detector}}) {
        my $detector_hash = $plan->{$detector}->{$version}->{'snv'}->[0];
        my $dir = "$build_dir/variants/snv/";
        $dir .= "$detector-" . $version . "-" . $detector_hash->{params};
        $dir =~ s/ /_/g;
        my $genome_build = $self->genome_build;
        Genome::Sys->shellcmd(cmd => "gmt vcf $module_name --output-file=$output_file --genome-build=$genome_build --type 'snv' --input-file $dir/snvs.hq --seq-center WUSTL --sample-id=$sample");
        if($self->annotate_dbsnp) {
            my $dbsnp_file = "$sample.$detector.dbsnp.vcf";
            Genome::Sys->shellcmd(cmd => "gmt vcf vcf-annotate-dbsnp --vcf-file=$output_file --output-file=$dbsnp_file --dbsnp_file = /gscuser/cmiller/snp130.txt");
            Genome::Sys->shellcmd(cmd => "mv $dbsnp_file $output_file");
        }
        my $filter_dir = "$dir/";
        for my $filter (@{$detector_hash->{filters}}) {
            $filter_dir .= $filter->{name} . "-" . $filter->{version} . "-" . $filter->{params} . "/";
            $filter_dir =~ s/ /_/g;
            my $filter_name = $filter->{name};
            my $filter_description = $filter->{class};
            my $filter_file = "$output_file.$filter_name";
            Genome::Sys->shellcmd(cmd => "gmt vcf vcf-filter --vcf-file=$output_file --output-file=$filter_file --filter-name=$filter_name --filter-description=$filter_description --filter-file=$filter_dir/snvs.hq --filter-keep --variant-type 'SNP'");
        Genome::Sys->shellcmd(cmd => "mv $filter_file $output_file");
        }
    }
        ###need to make a hash of result files
        $detector_vcf{$detector}=$output_file;
    }
    ###at this point we should have <number of detectors> vcf files laying around fully filtered and possibly dbsnp annotated. now merge them.
    my($merged_vcf, $complete_detector_list) =$self->depth_first_merge($tree->{'snv'}, $sample);
    my $output_file;

    ###now, apply LOH and DBSNP 
    if($build->loh_version) {
         $output_file = $self->output_dir . "/$sample.$complete_detector_list.merged.loh.vcf";
        Genome::Sys->shellcmd(cmd => "gmt vcf vcf-filter --vcf-file $merged_vcf --output-file $output_file --filter-name \"LOHfilter\" --filter-description \"LOH filter\" --filter-file $build_dir/loh/snvs.somatic.v2.bed --filter-keep --variant-type \"SNP\" --bed-input=1");
        Genome::Sys->shellcmd(cmd => "mv $output_file $merged_vcf");
    }
    if($build->previously_discovered_variations_build) {
         $output_file = $self->output_dir . "/$sample.$complete_detector_list.merged.dbsnp.vcf";
        Genome::Sys->shellcmd(cmd=> "gmt vcf vcf-filter --vcf-file $merged_vcf --output-file $output_file --filter-name \"novel\" --filter-description \"Novel Variant Filter\" --filter-file $build_dir/novel/snvs.hq.novel.v2.bed --filter-keep --variant-type \"SNP\" --bed-input=1");
        Genome::Sys->shellcmd(cmd => "mv $output_file $merged_vcf");
    }
    $dv2->delete();

    ###minimum work checking
    my %should_pass_hash;
    my @final_dv2_files = glob ("$build_dir/effects/snvs.hq.novel.tier[1234].v2.bed");
    unless(@final_dv2_files) {
        $self->status_message("Unable to QC my result-- cannot find $build_dir/effects/snvs.hq.novel.tier[1234].v2.bed");
        return 0;
    }
    for my $file (@final_dv2_files) {
        my $fh = IO::File->new($file);
        while(my $line = $fh->getline) {
            my ($chr, $start, $stop, $ref_var) = split /\t/, $line;
            $should_pass_hash{$chr}{$start+1}=1;
        }
    }
    my $fh = IO::File->new($merged_vcf);
    while(my $line = $fh->getline) {
        next if $line =~m/^#/;
        my ($chr, $pos, undef) = split /\t/, $line;
        if(exists($should_pass_hash{$chr}{$pos})) {
            if ($line =~ m/PASS/) {
                delete $should_pass_hash{$chr}{$pos};
            }
            else {
                $self->error_message("$line should be passing according to DV2 but is not");
            }
        }
    }
    for my $chr (keys %should_pass_hash) {
        for my $pos (keys %{$should_pass_hash{$chr}}) {
            $self->error_messagE("$chr, $pos expected in $merged_vcf according to DV2 but not found");
        }
    }



    return 1;
}


sub check_for_conversion_ability {
    my $self = shift;
    my $plan = shift;
    for my $detector (keys %{$plan}) {
        my $module = $detector_to_module{$detector};
        unless($module) {
            $self->error_message("Module to convert $detector not defined. Not creating vcf.");                  
            return 0;
        }
    }
    return 1;
}

sub depth_first_merge {
    my $self = shift;
    my $tree_part = shift;
    my $sample_name=shift;
    my @keys = keys %$tree_part;  
    #There should always be exactly one outer rule (or detector)
    unless(scalar @keys eq 1) {
        $self->error_message('Unexpected data structure encountered!  There were ' . scalar(@keys) . ' keys');
    }   

    my $key = $keys[0];

    #Case One:  We're somewhere in the middle of the data-structure--we need to combine some results
    if($key eq 'intersect' or $key eq 'union' or $key eq 'unionunique') {
        my $value = $tree_part->{$key};

        unless(ref $value eq 'ARRAY') {
            $self->error_message('Unexpected data structure encountered! I really wanted an ARRAY, not ' . (ref($value)||$value) );
            die($self->error_message);
        } 

        my ($keys_0, undef) = keys %{$value->[0]};
        my ($keys_1, undef) = keys %{$value->[1]};
        my ($file1, $file2, $merged_list1, $merged_list2);
        if($keys_0 ne 'detector') {
            ($file1, $merged_list1) =  $self->depth_first_merge($value->[0], $sample_name);
        }
        if($keys_1 ne 'detector') {
            ($file2, $merged_list2)  = $self->depth_first_merge($value->[1], $sample_name);
        }
        unless($file1) {
            $file1 = $detector_vcf{$value->[0]->{'detector'}->{'name'}};
        }
        unless($file2) {
            $file2 = $detector_vcf{$value->[1]->{'detector'}->{'name'}};
        } 
        my ($detector1, $detector2);
        if($merged_list1) {
            $detector1=$merged_list1;
        }  else {
            $detector1=$value->[0]->{'detector'}->{'name'};
        }
        if($merged_list2) {
            $detector2=$merged_list2;
        }       else {
            $detector2=$value->[1]->{'detector'}->{'name'};
        }

        

        my $output_file=$self->output_dir . "/$sample_name.$detector1.$detector2.merged.vcf";
        ###ok time to do an intersect, union, or union-unique
        $self->status_message("Merging $detector1 & $detector2-- $key");
        if($key eq 'intersect') {
            Genome::Sys->shellcmd(cmd => "gmt vcf vcf-merge --output-file=$output_file --vcf-files=$file1,$file2 --merge-filters --source-ids '$detector1, $detector2' --require-all-passing");
            return ($output_file, "$detector1.$detector2");
        }
        elsif($key eq 'union' || $key eq 'unionunique') {
            Genome::Sys->shellcmd(cmd=> "gmt vcf vcf-merge --output-file=$output_file --vcf-files=$file1,$file2 --merge-filters --source-ids '$detector1, $detector2' --keep-all-passing");
            return ($output_file, "$detector1.$detector2");
        }
    }
    else { #key not equal to intersect
        $self->error_message("expected an intersect or union operation at root of tree, got $key");
        return 0;
    }
}





1;

























