package Genome::Model::Tools::Analysis::ReportCoverage;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Sort::Naturally qw( nsort );

class Genome::Model::Tools::Analysis::ReportCoverage {
    is => 'Command',
    has => [
    build_ids => { 
        type => 'String',
        is_optional => 1,
        doc => "string of space separated build ids to check",
    },
    alignment_model_group => {
        type => 'String',
        is_optional => 1,
        doc => "The name of the Alignment Model Group to use",
    },
    somatic_model_group => {
        type => 'String',
        is_optional => 1,
        doc => "The name of the Somatic Model Group to use",
    },
    by_patient => {
        type => 'Boolean',
        is_optional => 1,
        doc => "Whether or not to report each patient as a single row instead of each genome",
        default => '0',
    },
    ]
};


sub execute {
    my $self=shift;
    $DB::single = 1;
    my @builds;
    if($self->build_ids) {
        my @build_ids = split / /, $self->build_ids;
        for my $build_id (@build_ids) {

            #grab the build
            my $build = Genome::Model::Build->get($build_id);
            unless(defined($build)) {
                $self->error_message("Unable to find build $build_id");
                return;
            }
            push @builds, $build;
        }
    }
    elsif($self->alignment_model_group) {
        my $group = Genome::ModelGroup->get(name => $self->alignment_model_group);
        unless($group) {
            $self->error_message("Unable to find a model group named " . $self->alignment_model_group);
            return;
        }
        for my $model ($group->models) {
            my $build = $model->last_complete_build;
            unless($build) {
                $self->error_message("No complete build for model " . $model->id);
            }
            else {
                push @builds, $build;
            }
        }
    }
    elsif($self->somatic_model_group) {
        my $group = Genome::ModelGroup->get(name => $self->somatic_model_group);
        unless($group) {
            $self->error_message("Unable to find a model group named " . $self->somatic_model_group);
            return;
        }
        for my $somatic_model ($group->models) {
            my $somatic_build = $somatic_model->last_succeeded_build or die "No succeeded build found for somatic model id $somatic_model.\n";
            my $somatic_build_id = $somatic_build->id or die "No build id found in somatic build object for somatic model id $somatic_model.\n";
            my $tumor_build = $somatic_build->tumor_build or die "Cannot find tumor build for somatic build $somatic_build_id.\n";
            my $normal_build = $somatic_build->normal_build or die "Cannot find normal build for somatic build $somatic_build_id.\n";
            push @builds, $tumor_build, $normal_build;
        }
    }
    else {
        $self->error_message("You must provide either build id(s) or a model group name to run this script");
        return;
    }
            
    my %coverage_metric_table;
    for my $build (@builds) {
        #grab metrics
        my %metrics = map {$_->name => $_->value} $build->metrics;

        #calculate common name like AML11
        my $common_name = $build->model->subject->source_common_name;

        #obtain build id
        my $build_id = $build->id;

        #this should tell us about whether it's tumor or normal
        my $type = $build->model->subject->common_name;
        my $subtype = $build->model->subject->sub_type;
        unless(defined $coverage_metric_table{$common_name}{subtype}) {
            $coverage_metric_table{$common_name}{subtype} = $subtype;
        }
        my @lanes = $build->instrument_data;
        $coverage_metric_table{$common_name}{$type} = {
            gbp => $metrics{"instrument data total kb"} ? sprintf("%0.01f",$metrics{"instrument data total kb"} / 10**6) : '-',
            haploid_coverage => $metrics{'haploid_coverage'} ? sprintf("%0.01f", $metrics{'haploid_coverage'}) : '-',
            unfiltered_diploid_coverage => $metrics{"gold-heterozygous-snp match heterozygous-1-allele-variant unfiltered"} ? sprintf("%0.02f",$metrics{"gold-heterozygous-snp match heterozygous-1-allele-variant unfiltered"}) : '-',
            filtered_diploid_coverage => $metrics{"gold-heterozygous-snp match heterozygous-1-allele-variant filtered"} ? sprintf("%0.02f",$metrics{"gold-heterozygous-snp match heterozygous-1-allele-variant filtered"}) : '-',
            lanes => scalar(@lanes),
            build_id => $build_id,
        };
            

    }

    if($self->by_patient) {
        #print header
        print "Sample\tSubtype\tTumor Lanes\tTumor Gbp\tTumor Haploid Coverage\tTumor Known Het SNP Coverage\tTumor Filtered Known Het SNP Coverage\tNormal Lanes\tNormal Gbp\tNormal Haploid Coverage\tNormal Known Het SNP Coverage\tNormal Filtered Known Het SNP Coverage\n";

        foreach my $sample (nsort  keys %coverage_metric_table) {
            my $subtype =  $coverage_metric_table{$sample}{subtype};
            print  "$sample\t$subtype";
            foreach my $type (qw( tumor relapse normal skin)) {
                next unless $coverage_metric_table{$sample}{$type};
                print  "\t", join("\t", @{$coverage_metric_table{$sample}{$type}}{qw( lanes gbp haploid_coverage unfiltered_diploid_coverage filtered_diploid_coverage)});
            }
            print "\n";
        }
    }
    else {
        print "Sample\tSubtype\tLanes\tGbp\tHaploid Coverage\tKnown Het SNP Coverage\tFiltered Known Het SNP Coverage\tBuild ID\n";

        foreach my $sample (nsort  keys %coverage_metric_table) {
            foreach my $type (qw( tumor relapse normal skin)) {
                next unless $coverage_metric_table{$sample}{$type};
                print  "$sample $type";
                print  "\t", join("\t", @{$coverage_metric_table{$sample}{$type}}{qw( lanes gbp haploid_coverage unfiltered_diploid_coverage filtered_diploid_coverage build_id)},), "\n";
            }
        }

    }

    return 1;

}


1;

sub help_brief {
    "Creates a table of coverages like you might frequently see at a monday meeting"
}

sub help_detail {
    <<'HELP';
HELP
}
