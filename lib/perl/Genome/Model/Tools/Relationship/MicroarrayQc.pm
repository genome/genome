package Genome::Model::Tools::Relationship::MicroarrayQc;

use strict;
use warnings;
use File::Basename;
use Genome;
use Genome::Info::IUB;
use Workflow;
use Workflow::Simple;

use List::MoreUtils qw/ uniq /;
class Genome::Model::Tools::Relationship::MicroarrayQc {
    is => 'Command',
    has => [
    output_dir => {
        is => 'Text',
        is_optional=>0,
        doc => 'the directory where you want results stored',
    },
    sample_list=> {
        is=>'Text',
        is_optional=>0,
        doc=>'TGI sample names file, one sample per line',
    },
    ref_build=> {
        is=>'Number',
        is_optional=>1,
        default=>106942997,
        doc=>"Do the comparison on build 37 iscan data by default"
    },
    ],
};

sub help_synopsis {
    return <<EOS
genome model-group relationship-qc --model-group=1745 --output-dir=/foo/bar/

EOS
}


sub execute {
    my $self=shift;
    $DB::single=1;
    my $output_dir = Genome::Sys->create_directory($self->output_dir);
    unless($output_dir) {
        $self->error_message("Unable to create output directory: " . $self->output_dir);
        return;
    }
    my $ref_build = Genome::Model::Build->get($self->ref_build);
    unless($ref_build) {
        $self->error_message("the supplied reference build id does not seem to be a real build.");
        return 0;
    }
    my $sample_fh= Genome::Sys->open_file_for_reading($self->sample_list);
    chomp(my @samples = $sample_fh->getlines);
    my %sample_to_model;
    my @missing_samples;
    ##gather all gold snp files for model
    for my $sample(@samples) {
        my @models = Genome::Model::GenotypeMicroarray->get(subject_name=>$sample, reference_sequence_build_id=>$self->ref_build);
        if(@models) {
            my $latest_model;
            for my $model (@models) {
                if($model->last_succeeded_build) {
                    if(!$latest_model || ($latest_model->id < $model->id)) {
                        $latest_model = $model;
                    }
                }
            }
            $sample_to_model{$sample}=$latest_model;
        }
        else {
            push @missing_samples, $sample;
        }
    }
    if(@missing_samples) {
        $self->status_message("#####MISSING MODELS:#########\n");
    }
    for my $missing (@missing_samples) {
        $self->error_message("No finished microarray models for $missing found on reference build " . $ref_build->name . "\n");
    }
    my $beagle_input = $self->generate_beagle_input_from_microarray(\%sample_to_model);  
    my $beagle_output = $self->run_beagle($beagle_input);
    $self->generate_relationship_table($beagle_output, $beagle_input);
    return 1;

}
sub generate_beagle_input_from_microarray {
    my $self = shift;
    my $sample_to_model = shift;
    my $out_file_name = $self->output_dir . "/iscan.bgl.input";
    my $output_file = Genome::Sys->open_file_for_writing($out_file_name);
    $self->status_message("#######FOUND MODELS:#######\n");
    my @header = qw| I ID|;
    my @fhs;
    for my $sample (sort keys %$sample_to_model) {
        my $model = $sample_to_model->{$sample};
        my $model_name = $model->name;
        my $model_id = $model->id;
        my $gold_snv_bed = $model->last_succeeded_build->snvs_bed;
        chomp(my $line_count = `wc -l $gold_snv_bed | cut -f1 -d' '`);
        $self->status_message("$sample: $model_id, $model_name: $line_count lines\n");
        push @fhs, Genome::Sys->open_file_for_reading($gold_snv_bed);
        push @header, ($sample, $sample);
    }
    my @current_lines = map {$_->getline} @fhs;
    ###write header
    $output_file->print(join("\t", @header) . "\n");
    while(my $id = find_matching_line(\@fhs, \@current_lines)) {
        my @output_line = ("M", $id);
        for my $line (@current_lines) {
            my ($chr, $start, $stop, $ref_var, undef) = split "\t", $line;
            my ($ref, $iub) = split "/", $ref_var;
            my $line_id = "$chr:$stop";
            my ($all1, $all2) = Genome::Info::IUB->iub_to_alleles($iub);
            unless($all1 && $all2) {
                die;
            }
            if($line_id ne $id) {
                die "lines are not synced: $line, $id\n";
            }
            push @output_line, ($all1, $all2);
        }
        $output_file->print(join("\t", @output_line) . "\n");
        @current_lines= map{$_->getline} @fhs;  
    }
    $output_file->close;
    return $out_file_name;
}









sub find_matching_line {
    my $fh_ref = shift;
    my $line_ref = shift;
    my $not_synced = 1;
    my @chrs;
    my @positions;
    for (my $i=0; $i < scalar(@{$line_ref}); $i++) {   
        my ($chr, $start, $pos, undef)  = split "\t", $line_ref->[$i];
        $chrs[$i]=$chr;
        $positions[$i]=$pos;
    }

    while(all_files_open($fh_ref)) {
        my $i = find_lowest_pos(@chrs);
        my $j = find_lowest_pos(@positions);
        if($i != -1) {
            $line_ref->[$i]=$fh_ref->[$i]->getline;
            my ($chr, $start, $pos, undef)  = split "\t", $line_ref->[$i];
            $chrs[$i]=$chr;
            $positions[$i]=$pos;
            next;
        }
        if($j !=-1) {
            $line_ref->[$j]=$fh_ref->[$j]->getline;
            my ($chr, $start, $pos, undef)  = split "\t", $line_ref->[$j];
            $chrs[$j]=$chr;
            $positions[$j]=$pos;
            next;
        }
        if(($i==-1) && ($j==-1)) {
            return $chrs[0] . ":".  $positions[0];
        }
    }
    return 0;
}







sub find_lowest_pos {
    my @positions = @_; 
    s/X/23/ for @positions;
    s/Y/24/ for @positions;
    my $least_position;
    my $least_idx;
    if(scalar(uniq @positions) ==1) {
        return -1; 
    }   
    for (my $i=0; $i < @positions; $i++) {
#        if($positions[$i]=~m/[MGL]/) {
#            return $i;
#        }
        if (!$least_position || $least_position > $positions[$i]) {
            $least_position = $positions[$i];
            $least_idx = $i; 
        }   
    }
    return $least_idx;
}



sub advance_a_line {
    my $fh_ref = shift;
    my $line_ref = shift;
    for (my $i = 0; $i < scalar(@{$fh_ref}); $i++) {
        $line_ref->[$i]=$fh_ref->[$i]->getline;
        if(!$line_ref->[$i]) {
            return 0;
        }
    }
    return 1;
}


sub all_files_open {
    my $fh = shift;
    for my $handle (@{$fh}) {
        if($handle->eof) {
            return 0;
        }
    }
    return 1;
}



sub run_beagle { 
    my $self = shift;
    my $beagle_input = shift;
    my $output_dir = $self->output_dir;

    #./beagle.sh  fastibd=true unphased=cleft_lip/mpileup/bgl_from_vcf.test_input out=cleft_lip/mpileup/cleft_lip.out missing=?
    my $cmd = "java -Xmx14000m -jar $ENV{GENOME_SW}/beagle/installed/beagle.jar fastibd=true unphased=$beagle_input out=$output_dir/beagle missing=?";
    my $rv = Genome::Sys->shellcmd(cmd=>$cmd, input_files=>[$beagle_input]);
    if($rv != 1) {
        $self->error_message("Error running Beagle\n");
        return;
    }
    my ($file, $dir) = fileparse($beagle_input);
    return $dir . "beagle." . $file . ".fibd.gz";
}

sub generate_relationship_table {
    my $self = shift;
    my $beagle_output = shift;
    my $fibd_file = Genome::Sys->open_gzip_file_for_reading($beagle_output);
    my $markers_file = shift;
    #TODO: make output files all in class def
    my $output_file = $self->output_dir . "/relationship_matrix.tsv";
    my $output_fh = Genome::Sys->open_file_for_writing($output_file);
    my $total_markers = `wc -l $markers_file`;
    $total_markers--; #account for header;
    my %relationships;
    while(my $line = $fibd_file->getline) {
        chomp($line);
        my ($first_guy, $second_guy, $start_marker, $stop_marker, $conf) = split "\t", $line;
        my $total_markers_covered;
        if(exists($relationships{$first_guy}{$second_guy})) {
            $total_markers_covered = $relationships{$first_guy}{$second_guy};
        }
        $total_markers_covered += ($stop_marker - $start_marker);
        $relationships{$first_guy}{$second_guy} = $total_markers_covered;
    }
    my $i=0;
    my %table_hash;
    my @table;
    my @header_row;
    $output_fh->print("INDIVIDUAL");
    for my $first_guy (sort keys %relationships) {
        $table_hash{$first_guy}=$i;
        $i++; 
        $output_fh->print("\t$first_guy");
        push @header_row, $first_guy;
    }
    for my $first_guy (sort keys %relationships) {
        for my $second_guy (sort keys %{$relationships{$first_guy}}) {
            unless(grep{/$second_guy/} @header_row) {
                $table_hash{$second_guy}=$i;
                $i++;
                $output_fh->print("\t$second_guy");
                push @header_row, $second_guy;
            }
        }
    }
    $output_fh->print("\n");



    for my $first_guy (sort keys %relationships) {
        for my $second_guy (sort keys %{$relationships{$first_guy}}) {
            my $total_shared_markers = $relationships{$first_guy}{$second_guy};
            my $percent = sprintf("%0.2f", $total_shared_markers/$total_markers * 100);
            my $j = $table_hash{$first_guy};
            my $k = $table_hash{$second_guy};
            $table[$j][$k]=$percent || 0;
        }
    }
    for (my $j=0; $j < $i; $j++) {
        my $person_for_row = $header_row[$j];
        $output_fh->print("$person_for_row");
        for (my $k=0; $k < $i; $k++) {
            my $value = $table[$j][$k] || "N/A";
            $output_fh->print("\t$value");
        }
        $output_fh->print("\n");
    }
    $output_fh->close();
    return 1;
}

1;
