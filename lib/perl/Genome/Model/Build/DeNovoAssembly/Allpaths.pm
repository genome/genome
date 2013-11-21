package Genome::Model::Build::DeNovoAssembly::Allpaths;;

use strict;
use warnings;
use Genome;
use List::Util;

class Genome::Model::Build::DeNovoAssembly::Allpaths {
    is => 'Genome::Model::Build::DeNovoAssembly',
};

sub validate_for_start_methods {
    my $self = shift;
    my @methods = $self->SUPER::validate_for_start_methods;
    push @methods, 'validate_sloptig_and_jump_instrument_data_assigned';
    return @methods;
}

sub validate_sloptig_and_jump_instrument_data_assigned {
    my $self = shift;

    my ($jumping_count, $sloptig_count) = (qw/ 0 0 /);
    for my $i_d ($self->instrument_data) {
        if ($self->_instrument_data_is_jumping($i_d)) {
            $jumping_count++;
        }
        elsif ($self->_instrument_data_is_sloptig($i_d)) { 
            $sloptig_count++;
        }
    }

    if ($jumping_count == 0) {
        return UR::Object::Tag->create(
            properties => ['instrument_data'],
            desc => "No jumping library instrument data found",
        );
    }

    if ($sloptig_count == 0) {
        return UR::Object::Tag->create(
            properties => ['instrument_data'],
            desc => "No sloptig library instrument data found",
        );
    }
    
    return;
}

sub _instrument_data_is_jumping {
    my ($self, $instrument_data) = @_;
    if ($self->_instrument_data_is_valid($instrument_data) and !($self->_instrument_data_is_sloptig($instrument_data))) {
        return 1;
    }
    return;
}

sub _instrument_data_is_sloptig {
    my ($self, $instrument_data) = @_;
    if ($self->_instrument_data_is_valid($instrument_data) and $instrument_data->original_est_fragment_size <=
        2.5*$instrument_data->read_length){
        return 1;
    }
    return;
}

sub _sx_result_is_x {
    my ($self, $sx_result, $x) = @_;
    my $answer;
    my $count = 0;
    my $func_name = "_instrument_data_is_".$x;
    for my $id ($sx_result->instrument_data) {
        my $temp = $self->$func_name($id);
        if ($count == 0) {
            $answer = $temp;
        }
        elsif ($answer == $temp) {
            next;
        }
        else {
            $self->error_message("Some instrument data in the sx result were sloptigs and some were not");
            die $self->error_message;
        }
    }
    return $answer;
}

sub _instrument_data_is_valid {
    my ($self, $instrument_data) = @_;
    unless($instrument_data->original_est_fragment_size and $instrument_data->read_length and
    $instrument_data->original_est_fragment_std_dev and $instrument_data->read_orientation) {
        my $error = join(" ","Instrument data did not have all required fields set: ",
                         "id:",
                         $instrument_data->id,
                         "original_est_fragment_size:",
                         $instrument_data->original_est_fragment_size,
                         "read_length:",
                         $instrument_data->read_length,
                         "original_est_fragment_std_dev:",
                         $instrument_data->original_est_fragment_std_dev,
                         "read_orientation: ",
                         $instrument_data->read_orientation
                         );
                         
        $self->error_message($error);
        return;
    }
    return 1;
}

#Override base class method
sub stats_file {
    my $self = shift;
    return $self->data_directory."/metrics.out";
}

sub _allpaths_in_group_file {
    return $_[0]->data_directory."/in_group.csv";
}

sub _allpaths_in_libs_file {
    return $_[0]->data_directory."/in_libs.csv";
}

sub before_assemble {
    my $self = shift;
    my @sx_results = @_;
    $self->status_message("Allpaths config files");

    my %params = $self->processing_profile->assembler_params_as_hash;

    $self->status_message("Generating Allpaths in_group.csv and in_libs.csv");
    my $in_group = "file_name,\tlibrary_name,\tgroup_name";

    my $in_libs = "library_name,\tproject_name,\torganism_name,\ttype,\tpaired,\tfrag_size,\tfrag_stddev,\tinsert_size,\tinsert_stddev,\tread_orientation,\tgenomic_start,\tgenomic_end";

    my %libs_seen;
    for my $sx_result (@sx_results) {
        my @lib_ids = $self->resolve_attribute_for_instrument_data("library_id", 0, $sx_result->instrument_data);
        my $lib_id = join("-", @lib_ids);
        my $genomic_start = 0;
        my $genomic_end = 0;
        my $lib_name = join("-", map{Genome::Library->get($_)->name} @lib_ids);
        if ($lib_name =~ /CHORI/) {
            $genomic_start = 2;
            $genomic_end = 50;
        }

        $in_group = $in_group."\n".$self->data_directory."/".$sx_result->resolve_base_name_from_instrument_data($sx_result->instrument_data).".*.".$sx_result->output_file_suffix.",\t".$lib_name.",\t".$sx_result->id;
        
        if (! $libs_seen{$lib_id}) {
            my $orientation;
            my $read_orientation = $self->resolve_attribute_for_instrument_data("read_orientation", 1, $sx_result->instrument_data);
            if ($read_orientation eq "forward_reverse") {
                $orientation = "inward";
            }
            elsif ($read_orientation eq "reverse_forward") {
                $orientation = "outward";
            }
            else {
                my @instrument_data_id_list = $self->resolve_attribute_for_instrument_data("id", 0, $sx_result->instrument_data);
                my $instrument_data_ids = join(",", @instrument_data_id_list);
                my $error = join(" ","Instrument data with ids", $instrument_data_ids, "has unrecognized read orientation", $read_orientation);
                $self->error_message($error);
                return;
            }
            
            my $species_name = $self->resolve_attribute_for_instrument_data("species_name", 1, $sx_result->instrument_data);
            my $original_est_fragment_std_dev = $self->resolve_average_for_attribute(attribute => "original_est_fragment_std_dev", objects => [$sx_result->instrument_data]);
            my $original_est_fragment_size = $self->resolve_average_for_attribute(attribute => "original_est_fragment_size", objects => [$sx_result->instrument_data]);
            if ($self->_sx_result_is_x($sx_result, "sloptig")) {
                my $fragment_std_dev = $original_est_fragment_std_dev;
                $in_libs = $in_libs."\n".$lib_name.",\tproject_name,\t".$species_name.",\tfragment,\t1,\t".$original_est_fragment_size.",\t".$fragment_std_dev.",\t,\t,\t".$orientation.",\t$genomic_start,\t$genomic_end";
            }
            elsif ($self->_sx_result_is_x($sx_result, "jumping")){
                my $fragment_std_dev = $original_est_fragment_std_dev;
                $in_libs = $in_libs."\n".$lib_name.",\tproject_name,\t".$species_name.",\tjumping,\t1,\t,\t,\t".$original_est_fragment_size.",\t".$fragment_std_dev.",\t".$orientation.",\t$genomic_start,\t$genomic_end";
            }
        }
        $libs_seen{$lib_id} = 1;
    }

    my $in_group_file = $self->_allpaths_in_group_file;
    unlink $in_group_file if -e $in_group_file;
    $self->status_message("Allpaths in_group file: ".$in_group_file);
    my $fh = eval { Genome::Sys->open_file_for_writing( $in_group_file); };
    if (not $fh) {
        $self->error_message("Can not open file ($in_group_file) for writing $@");
        return;
    }
    $fh->print($in_group);
    $fh->close;
    $self->status_message("Allpaths in_group file...OK");

    my $in_libs_file = $self->_allpaths_in_libs_file;
    unlink $in_libs_file if -e $in_libs_file;
    $self->status_message("Allpaths in_libs file: ".$in_libs_file);
    $fh = eval { Genome::Sys->open_file_for_writing( $in_libs_file); };
    if (not $fh) {
        $self->error_message("Can not open file ($in_libs_file) for writing $@");
        return;
    }
    $fh->print($in_libs);
    $fh->close;
    $self->status_message("Allpaths in_libs file...OK");
}

sub assembler_params {
    my $self = shift;

    my %default_params = (
        run => "run",
        sub_dir => "test",
        reference_name => "sample",
    );
    my %params = $self->processing_profile->assembler_params_as_hash;

    for my $param (keys %default_params) {
        if (! defined $params{$param}) {
            $params{$param} = $default_params{$param};
        }
    }

    $params{version} = $self->processing_profile->assembler_version;
    $params{pre} = $self->data_directory;
    $params{in_group_file} = $self->_allpaths_in_group_file;
    $params{in_libs_file} = $self->_allpaths_in_libs_file;
    $params{max_memory_gb} = $self->_mem_in_gb;

    return %params;
}

sub read_processor_output_file_count_for_instrument_data {
    my $self = shift;
    my $instrument_data = shift;

    if ($instrument_data->is_paired_end) {
        return 2;
    }
    else {
        return 1;
    }
}

sub read_processor_params_for_instrument_data {
    my $self = shift;
    my $instrument_data = shift;

    my $read_processor = $self->processing_profile->read_processor;

    my $output_file_count = $self->read_processor_output_file_count_for_instrument_data(    $instrument_data);

    return (
        instrument_data_id => $instrument_data->id,
        read_processor => $read_processor,
        output_file_count => $output_file_count,
        output_file_type => 'sanger',
        test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
    );
}

sub resolve_assemble_lsf_resource {
    my $self = shift;

    my $mem = $self->_mem_in_gb;

    my $mem_reserve = $mem*1024;
    my $mem_limit = $mem_reserve*1024 + 1048576;
    my $template = "-n 4 -R 'span[hosts=1] select[type==LINUX64 && mem>%s] rusage[mem=%s]' -M %s";
    return sprintf($template, $mem_reserve, $mem_reserve, $mem_limit);
}

sub _mem_in_gb {
 
    my $self = shift;
    my $mem = 494;
    my $egs = $self->_get_estimated_genome_size();

    if ($ENV{UR_DBI_NO_COMMIT}) {
        $mem = 60;
    } elsif ($egs and $egs <= 40_000_000) {
        $mem = 200;
    }
    return $mem;
}

sub resolve_assemble_lsf_queue {
    my $self = shift;
    my $queue = 'assembly';
    $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_PROD} if $ENV{UR_DBI_NO_COMMIT};
    return $queue;
}

sub calculate_estimated_kb_usage {
    my $self = shift;
    my $kb_per_read_count = 2; #based on the first 15 succeeded builds

    my $min_kb_reserved = 400_000_000;
    my $total_read_count = 0;
    for my $id ($self->instrument_data) {
        $total_read_count += $id->read_count;
    }
    my $estimate = $total_read_count*$kb_per_read_count;

    return List::Util::max($min_kb_reserved, $estimate);
}

sub _get_estimated_genome_size {
    my $self = shift;
    if ($self->can('subject')) {
        my $subject = $self->subject;
        return unless ($subject);
        if ($subject->can('estimated_genome_size')) {
            return $subject->estimated_genome_size;
        }
        if ($subject->can('taxon')) {
            my $taxon = $subject->taxon;
            return unless $taxon;
            return $taxon->estimated_genome_size;
        }
    }
    return;
}

1;
