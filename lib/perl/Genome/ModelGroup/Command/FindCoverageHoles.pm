package Genome::ModelGroup::Command::FindCoverageHoles;

use strict;
use warnings;

use Genome;
use Workflow;
use Workflow::Simple;
use List::Util qw(sum);
class Genome::ModelGroup::Command::FindCoverageHoles {
    is => 'Genome::Command::Base',
    has => [
    model_group_id=> {
        is=>'String',
        is_input=>1,
        is_optional=>0,
    },
    list_of_regions_to_examine=> { 
        is=>'String',
        is_input=>1,
        is_optional=>0,
    },
    coverage_cutoff=> {
        is=>'Number',
        is_input=>1,
        is_optional=>1,
        default=>5,
    },
    annotation_substructures_bed=> {
        is=>'String',
        default=>"/gscmnt/gc4095/info/feature_list/A8375E4058CD11E196940D538BF3D0F7/A8375E4058CD11E196940D538BF3D0F7.bed",
        doc=>"this is for Ensembl only build 37. it should be gettable for any build w/ Genome::Model::Build->get(id=>106409296)->get_or_create_roi_bed->file_path, but there's a bug in that that makes new allocations each time is run",
    },
    output_dir=> {
        is=>'String',
        is_input=>1,
        is_output=>1,
        is_optional=>0,
        doc=>'the output directory',
    },
    cleanup=> {
        is=>'Number',
        default=>1,
        doc=>"if set, will leave behind a 'by_bases' dir with all the readcounts. you could use this to graph coverage pretty easily'"
    },
    ],
    has_optional => {
    },
    doc => 'make a new model group from another, varying properties on the model as specified' 
};

sub help_synopsis {
    return <<EOS
genome model-group find-coverage-holes --model-group-id=80085 --list-of-regions=some_regions.bed

EOS
}


sub execute {
    my $self = shift;
    $DB::single=1;
    my $model_group_id = $self->model_group_id;
    $self->status_message("Model Group Id set to: $model_group_id...getting...");
    my $mg = Genome::ModelGroup->get($model_group_id);
    unless($mg) {
        $self->error_message("Unable to ->get a model group with id $model_group_id");
        return 0;
    }
    my %bams = map { $_->subject->name =>$_->last_succeeded_build->whole_rmdup_bam_file } $mg->models;
    my $input_bed = $self->list_of_regions_to_examine;
    my $bed = $self->copy_bed_and_add_region_tags($self->output_dir, $input_bed);
    unless(-d $self->output_dir) {
        Genome::Sys->create_directory($self->output_dir);
    } 
     my @readcount_files = $self->generate_by_base_coverage($self->output_dir, \%bams, $bed);
    my $cutoff = $self->coverage_cutoff;
    my $bad_regions_bed = $self->find_low_coverage_regions($self->output_dir, $cutoff, \@readcount_files);
    my $intersected_file = $self->intersect_annotation_with_bad_regions($self->output_dir, $bad_regions_bed, $self->annotation_substructures_bed);
    $self->trim_and_pretty_output($self->output_dir, $intersected_file);
    #cleanup
    if($self->cleanup) {
        unlink(@readcount_files);
        rmdir($self->output_dir . "/by_base/");
    }
    return 1;
}

sub copy_bed_and_add_region_tags {
    my ($self, $output_dir, $input_bed) = @_;
    my $output_bed = "$output_dir/bed_used_for_analysis.bed";
    my $fh = IO::File->new($input_bed);
    my $line = $fh->getline;
    chomp($line);
    my ($chr, $start, $stop, $region) = split "\t", $line;
    
    if(!defined($region)) {
        my $region_count=0;
        my $ofh = IO::File->new($output_bed, ">");
        $ofh->print($line . "\tr$region_count\n");
        while(my $line = $fh->getline) {
            $region_count++;
            chomp($line);
            $ofh->print("$line\tr$region_count\n");
        }
        $ofh->close();
    }else {
        `cp $input_bed $output_bed`;
    }
    return $output_bed;
}




sub trim_and_pretty_output {
    my ($self, $output_dir, $intersected_file) = @_;
    my $fh = Genome::Sys->open_file_for_reading($intersected_file);
    my %uncovered;
    my $output_file = $output_dir . "/annotated_poorly_covered_regions.bed";
    while(my $line = $fh->getline) {
        chomp($line);
        my ($chr, $start, $stop, $coverage, $other_chr, $other_start, $other_pos, $gene_transcript, $overlap) = split "\t", $line;
        my ($gene, undef) = split ":", $gene_transcript;
        if(exists($uncovered{$chr}{$start})) { 
            my $other_gene = $uncovered{$chr}{$start}{gene};
            my @genes = split ":", $other_gene;
            if(!(grep {/$gene/} @genes)) {
                my $current_line =  $uncovered{$chr}{$start}{line};
                my ($cchr, $cstart, $cstop, $ccoverage, $cgene_transcript) = split "\t", $current_line;
                $cgene_transcript.="/$gene_transcript";
                $current_line = "$cchr\t$cstart\t$cstop\t$ccoverage\t$cgene_transcript";     
                $uncovered{$chr}{$start}{gene}.="$gene:";
                $uncovered{$chr}{$start}{line}=$current_line;
            }
        }else{
            $uncovered{$chr}{$start}{gene}="$gene:";                       
            my $line_to_store = "$chr\t$start\t$stop\t$coverage\t$gene_transcript"; 
            $uncovered{$chr}{$start}{line}=$line_to_store;
        }
    }
    my $ofh = IO::File->new($output_file, ">");
    for my $chr (sort keys %uncovered) {
        for my $start (sort keys %{$uncovered{$chr}}) {
            my $line = $uncovered{$chr}{$start}{line};
            $ofh->print("$line\n");
        }
    }
}

sub intersect_annotation_with_bad_regions {
    my ($self, $output_dir, $bad_regions_bed, $annotation_bed) = @_;
    my $temp = Genome::Sys->create_temp_file_path();
    my $cmd = "intersectBed -a $bad_regions_bed -b $annotation_bed -wao > $temp";
    Genome::Sys->shellcmd(cmd=>$cmd);
    return $temp;
}

sub open_gzip_files {
    my ($self, $readcount_files) = @_;
    my @opened_file_handles;
    for my $file (@$readcount_files) {
        my $fh = Genome::Sys->open_gzip_file_for_reading($file);
        push @opened_file_handles, $fh;
    }
    return \@opened_file_handles;
}


sub find_low_coverage_regions {
    my ($self, $output_dir, $cov_cutoff, $readcount_files) = @_;
    my $output_file = $output_dir . "/poorly_covered_regions.bed";
    my $out_fh = IO::File->new($output_file, ">");

    my $by_base_files = $self->open_gzip_files($readcount_files);

    my($start_chr, $start_pos, $start_region_id,$total_depth)=(0,0,0,0);
    my($last_chr, $last_pos, $last_region_id);
    my @last_depths;
    my ($chr, $pos, $current_region_id, $average_depth) =$self->process_readcount_lines($by_base_files);
    while($chr) {
        if($current_region_id ne $last_region_id) {
            @last_depths=();
        }
        push @last_depths, $average_depth;
        if(scalar(@last_depths) > 30) {
            shift @last_depths;
        }
        my $windowed_average  = sum(@last_depths) / scalar(@last_depths);
        if($windowed_average <= $cov_cutoff) {

            $total_depth+=$average_depth;
            if($start_pos) {
                if($current_region_id ne $start_region_id) {
                    #stop current region and write out, then start a new one
                    $self->write_region($out_fh, $start_chr, $start_pos, $last_pos, $total_depth);
                    ($start_chr, $start_pos, $start_region_id,$total_depth)=($chr, $pos, $current_region_id, $average_depth);
                    @last_depths=();
                }
                else {
                    #do nothing, region will be extended
                }
            }
            else { #no region has been started, start one
                ($start_chr, $start_pos, $start_region_id, $total_depth) = ($chr, $pos, $current_region_id, $average_depth);
            }
        } else { #position is greater than cutoff, so regions should end, or nothing should happen
            if($start_pos) {
                if($current_region_id ne $start_region_id) {
                    #stop current region and write out
                    $self->write_region($out_fh, $start_chr, $start_pos, $last_pos, $total_depth);
                    @last_depths=();
                }
                else {
                    $self->write_region($out_fh, $start_chr, $start_pos, $pos, $total_depth);
                }
                ($start_chr, $start_pos, $start_region_id, $total_depth) = (0,0,0,0);
            }
        }
        ($last_chr, $last_pos, $last_region_id) = ($chr, $pos, $current_region_id);
        ($chr, $pos, $current_region_id, $average_depth) =$self->process_readcount_lines($by_base_files);
    }
    return $output_file;
}

sub process_readcount_lines {
    my ($self, $files) = @_;
    my ($return_chr, $return_pos, $current_region_id, $average_depth);
    for my $file (@$files) {
        my $line = $file->getline;
        if(!$line) {
            return undef;
        }
        chomp($line);
        #expecting: 1   18772299    19208054    r0  1   83
        my ($chr, $reg_start, $reg_stop, $region_id, $offset, $coverage) = split "\t", $line;
        #might have gotten 1   18772299    19208054    1   83 if the region names weren't in input bed.
        unless(defined($coverage)) {
            $coverage=$offset;
            $offset=$region_id;
        }


        if($return_pos) {
            die unless ($reg_start+$offset)==($return_pos);
        }
        else {
            $return_chr = $chr;
            $return_pos = $reg_start+$offset;
            $current_region_id=$region_id;
        }
        $average_depth+=$coverage;
    }
    if($average_depth > 0) {
        $average_depth/=scalar(@$files);
    } 
    return ($return_chr, $return_pos, $current_region_id, $average_depth);
}








sub write_region {
    my ($self, $out_fh, $start_chr, $start_pos, $stop_pos, $total_depth)= @_;
    my ($chr, $start, $stop, $region_depth);
    $chr = $start_chr;
    $start = $start_pos -1;
    $stop = $stop_pos;
    $region_depth = ($total_depth) / ($stop - $start);
    $out_fh->print("$chr\t$start\t$stop\t$region_depth\n");
}





sub generate_by_base_coverage {
    my ($self, $output_dir, $bams, $bed) = @_;
    my $by_base_dir = $output_dir . "/by_base/";
    unless(-d $by_base_dir) {
        Genome::Sys->create_directory($by_base_dir);
    }


    my %input_properties;
    my @inputs=("bed_file");
    my @output_files;
    $input_properties{bed_file}=$bed;

    for my $sample (sort keys %$bams) {
        my $bam = $bams->{$sample};
        $input_properties{"bam_$sample"}=$bam;
        my $output_file = $by_base_dir . "$sample.bybase.output";
        push @output_files, $output_file;
        $input_properties{"output_$sample"}=$output_file;
        push @inputs, ("bam_$sample", "output_$sample");
    }
    my $number = scalar(keys %$bams);
    my $workflow = Workflow::Model->create(
        name=> "$number sample parallel coverage generation",
        input_properties => [
        @inputs,
        ],
        output_properties => [
        'output',
        ],
    );
    $workflow->log_dir($self->output_dir);
    for my $sample (sort keys %$bams) {
        my $coveragebed_op = $workflow->add_operation(
            name=>"coverageBed $sample",
            operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Bed::CoverageBed"),
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"bed_file",
            right_operation=>$coveragebed_op,
            right_property=>"bed_file",
        );

        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"bam_$sample",
            right_operation=>$coveragebed_op,
            right_property=>"bam_file",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"output_$sample",
            right_operation=>$coveragebed_op,
            right_property=>"output_file",
        );
        $workflow->add_link(
            left_operation=>$coveragebed_op,
            left_property=>"output_file",
            right_operation=>$workflow->get_output_connector,
            right_property=>"output",
        );
    }
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }
    $self->status_message("Now launching $number coverageBed jobs");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %input_properties);
    unless($result) {
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("parallel mpileup workflow did not return correctly.");
    }
    return @output_files;


}









1;
