package Genome::SoftwareResult::Command::FixLookupHashesDV2;

use strict;
use warnings;

use Genome;
use Genome::Utility::List "in";
use JSON;


class Genome::SoftwareResult::Command::FixLookupHashesDV2 {
    is => 'Command::V2',

    has_optional => [
        commit_size => {
            is => 'Number',
            default_value => '100',
            doc => 'Commit to the UR::Context with this many ' .
                'calculated hashes at a time.',
        },
        class_name => {
            is => 'Text',
            valid_values => [
                'Genome::Model::Tools::DetectVariants2::Result::Vcf::Combine',
                'Genome::Model::Tools::DetectVariants2::Result::Vcf::Filter'],
            doc => 'Fill or replace the lookup_hash values for ' .
                'all SoftwareResults of this class.',
        },
        num_cpus => {
            is => 'Number',
            default_value => 1,
            doc => 'Use this many processors to fill.',
        },
        output_directory => {
            is => 'Text',
        },
        blacklist_id_filename => {
            is => 'Text',
            default_value => 'blacklist_id_file.txt',
        },
        debug => {
            is => 'Boolean',
            default_value => 0,
        },
    ],

};

sub execute {
    my $self = shift;

    my $previous_progress_index = $self->_read_progress_index('main');
    my @pids;
    for (my $i = 0; $i < $self->num_cpus; $i++) {
        my $pid = UR::Context::Process->fork();
        unless (defined($pid)) {
            die $self->error_message("Failed to fork on subprocess number $i");
        }

        if ($pid == 0) {
            # child
            $self->fill_hashes($self->num_cpus, $i);
            exit;
        } else {
            # parent
            push(@pids, $pid);
        }
    }

    for my $pid (@pids) {
        waitpid($pid, 0);
    }

    $self->_combine_all_progress_indecies($self->num_cpus);
    return 1;
}

sub _is_combine_type {
    my ($self) = @_;
    if ($self->class_name =~ m/Combine/) {
        return 1;
    } else {
        return 0;
    }
}

sub _progress_index_filename {
    my ($self, $remainder) = @_;

    my $type;
    if ($self->_is_combine_type) {
        $type = 'combine';
    } else {
        $type = 'filter';
    }
    my $filename = $type . "_$remainder" . "_progress_index.txt";
    return File::Spec->join($self->output_directory, $filename);
}

sub _read_progress_index {
    my ($self, $remainder) = @_;

    my $filename = $self->_progress_index_filename($remainder);
    my $progress_index;
    if (-e $filename) {
        $self->status_message("Reading in $filename.");
        my $json = new JSON;
        my $fh = Genome::Sys->open_file_for_reading($filename);
        my $json_text = <$fh>;
        $progress_index = $json->decode($json_text);
        $self->status_message(
            sprintf("Read in %d filled_ids, and %d unfillable_ids",
            scalar(keys %{$progress_index->{'filled_ids'}}),
            scalar(keys %{$progress_index->{'unfillable_ids'}})));
    } else {
        $self->status_message("Couldn't find $filename to open as json text.");
        $progress_index = {filled_ids => {}, unfillable_ids => {}};
    }
    # debug how many read in and such.
    return $progress_index;
}

sub _write_progress_index {
    my ($self, $progress_index, $remainder) = @_;

    my $json = new JSON;
    my $json_text = $json->encode($progress_index);

    my $filename = $self->_progress_index_filename($remainder);
    my $fh = IO::File->new($filename, 'w');
    print $fh $json_text;
    $fh->close();
}

sub _combine_all_progress_indecies {
    my ($self, $base) = @_;

    my $main_progress_index = $self->_read_progress_index('main');
    for my $i (0..$base-1) {
        my $progress_index = $self->_read_progress_index($i);
        for my $key (keys %$progress_index) {
            my $main_one = $main_progress_index->{$key};
            my $this_one = $progress_index->{$key};
            for my $sub_key (keys %$this_one) {
                $main_one->{$sub_key} = $this_one->{$sub_key};
            }
        }
    }
    $self->_write_progress_index($main_progress_index, 'main');
}

sub _should_skip {
    my ($self, $base, $remainder, $progress_index, $sr_id) = @_;

    return 1 unless $sr_id % $base == $remainder;
    return 1 if exists $progress_index->{filled_ids}->{$sr_id};
    return 1 if exists $progress_index->{unfillable_ids}->{$sr_id};
    return 0;
}

sub fill_hashes {
    my ($self, $base, $remainder, $previous_progress_index) = @_;

    my %progress_index;
    my $total = 0;
    my $i = 0;
    my $iterator = $self->_get_iterator();
    my $progress_index = $self->_read_progress_index($remainder);

    while (my $sr = $iterator->next()) {
        my $sr_id = $sr->id;
        next if $self->_should_skip($base, $remainder, $previous_progress_index, $sr_id)
                or $self->_should_skip($base, $remainder, $progress_index, $sr_id);
        #$self->status_message("Not skipping $sr_id.");

        my $fixed_status = 0;
        if ($self->_is_combine_type) {
            $fixed_status = fix_combine_result($sr, $self->debug);
        } else {
            $fixed_status = fix_filter_result($sr, $self->debug);
        }
        if ($fixed_status) {
            $sr->lookup_hash($sr->calculate_lookup_hash());
            $progress_index->{'filled_ids'}->{$sr_id} = 1;
            $i++;
        } else {
            $progress_index->{'unfillable_ids'}->{$sr_id} = 1;
        }

        if ($i >= $self->commit_size) {
            $total += $i;
            $self->status_message("Committing $i of $total rows on process $remainder");
            $i = 0;

            # Destroy iterator before commit, so we don't get a warning
            undef $iterator;
            UR::Context->commit();

            $self->_write_progress_index($progress_index, $remainder);

            $iterator = $self->_get_iterator();
        }
    }

    if ($i) {
        $total += $i;
        $self->status_message("Committing $i of $total rows on process $remainder");

        # Destroy iterator before commit, so we don't get a warning
        undef $iterator;
        UR::Context->commit();
        $self->_write_progress_index($progress_index, $remainder);
    }
}

sub _get_iterator {
    my ($self) = @_;

    my $iterator;
    if ($self->class_name) {
        $iterator = Genome::SoftwareResult->create_iterator(
            'subclass_name' => $self->class_name);
    }
    return $iterator;
}


sub fix_filter_result {
    my ($result_to_fix, $debug) = @_;

    print "Working with result id " . $result_to_fix->id . "\n" if $debug;

    # Get the input result, the one that this result translated to vcf
    my $pre_vcf_conversion_input = $result_to_fix->input;
    die "Could not get input result for id " . $result_to_fix->id . " and input id " . $result_to_fix->input_id . "\n" unless $pre_vcf_conversion_input;

    # Now that I have the input to this result, get ITS input to find the step before this one
    # FIXME is this at all correct?

    # Get the things that I am using (the result that came before me)
    my @user_bridges = Genome::SoftwareResult::User->get(user => $pre_vcf_conversion_input);
    return 0 unless @user_bridges;
    print "Too many user bridges: " . join(',', map {$_->id} @user_bridges) . "\n" if (@user_bridges > 1);
    return 0 if (@user_bridges > 1);

    my $result_i_am_filtering = $user_bridges[0]->software_result;

    # Get the vcf result for the non-vcf step before this. This will be the vcf file that came before this vcf file.
    # ...with OUR vcf version (to connect each version together as appropriate)
    my $aligned_reads_sample = $result_to_fix->aligned_reads_sample;
    my $vcf_result_i_am_filtering = get_vcf_result($result_i_am_filtering, $result_to_fix->vcf_version,
            $aligned_reads_sample) || return 0;

    # Set the current result from the non-vcf result
    my $new_name = $pre_vcf_conversion_input->filter_name;
    my $new_params = $pre_vcf_conversion_input->filter_params || " "; # FIXME ... why does blank string not work?
    my $new_version = $pre_vcf_conversion_input->filter_version;
    my $new_strat = $pre_vcf_conversion_input->previous_filter_strategy || " ";
    if ($debug) {
        print "Non vcf input to this step is " . $pre_vcf_conversion_input->id . "\n";
        print "Setting incoming vcf result to " . $vcf_result_i_am_filtering->id . " ... the vcf result of the non vcf input\n";
        print "Setting filter name to $new_name\n";
        print "Setting filter params to $new_params\n";
        print "Setting filter version to $new_version\n";
        print "Setting previous filter strategy to $new_strat\n";
    }
    # Set the current result (which is a vcf result)'s incoming_vcf_result to the thing we found
    $result_to_fix->incoming_vcf_result($vcf_result_i_am_filtering);

    $result_to_fix->filter_name($new_name);
    $result_to_fix->filter_params($new_params);
    $result_to_fix->filter_version($new_version);
    $result_to_fix->previous_filter_strategy($new_strat);

    my @users = _get_user_objects($vcf_result_i_am_filtering);
    unless (in($result_to_fix,@users)) {
        print "Sanity check failed... I am not using the previous result's vcf result " .
            $result_to_fix->id . "\n";
        return 0;
    }
    return 1;
}

sub _get_user_objects {
    my ($sr) = @_;

    my @users = map{$_->user} $sr->users;
    return @users;
}

sub fix_combine_result {
    my ($result_to_fix, $debug) = @_;

    print "Working with result id " . $result_to_fix->id . "\n" if $debug;
    # Get the input result, the one that this result translated to vcf
    my $pre_vcf_conversion_input = $result_to_fix->input;
    die "Could not get input result for id " . $result_to_fix->id . " and input id " . $result_to_fix->input_id . "\n" unless $pre_vcf_conversion_input;

    # Now that I have the input to this result, get ITS inputs to find the steps before this one
    my $result_i_am_combining_a = $pre_vcf_conversion_input->input_a || die "cannot find input_a\n";
    my $result_i_am_combining_b = $pre_vcf_conversion_input->input_b || die "cannot find input_b\n";

    # get the vcf result of the previous result with OUR vcf version (to connect each version together as appropriate)
    my $input_result_vcf_a = get_vcf_result($result_i_am_combining_a, $result_to_fix->vcf_version) || return 0;
    my $input_result_vcf_b = get_vcf_result($result_i_am_combining_b, $result_to_fix->vcf_version) || return 0;

    my $input_a_id = $result_i_am_combining_a->id;
    my $input_b_id = $result_i_am_combining_b->id;

    if ($debug) {
        print "The (non-vcf) input for this result is " . $pre_vcf_conversion_input->id . "\n";
        print "Input a (non vcf) for the non-vcf result is " . $result_i_am_combining_a->id . "\n";
        print "Input b (non vcf) for the non-vcf result is " . $result_i_am_combining_b->id . "\n";
        print "The id for input a is $input_a_id... setting input_a_id to this\n";
        print "The id for input b is $input_b_id... setting input_b_id to this\n";
        print "The vcf result id for input a is " . $input_result_vcf_a->id . " ... setting incoming vcf result a to this\n";
        print "The vcf result id for input b is " . $input_result_vcf_b->id . " ... setting incoming vcf result b to this\n";
    }

    my @users_a = _get_user_objects($input_result_vcf_a);
    my @users_b = _get_user_objects($input_result_vcf_b);
    unless (in($result_to_fix,@users_a) and in($result_to_fix,@users_b)) {
        print "Sanity check failed... I am not using the previous result's vcf result " .
            $result_to_fix->id . "\n";
        return 0;
    }

    my @nonvcf_users_a = _get_user_objects($result_i_am_combining_a);
    my @nonvcf_users_b = _get_user_objects($result_i_am_combining_b);
    unless (in($pre_vcf_conversion_input,@nonvcf_users_a) and in($pre_vcf_conversion_input,@nonvcf_users_b)) {
        print "Sanity check failed... my nonvcf result is not using the previous nonvcf results" .
            $result_to_fix->id . "\n";
        return 0;
    }

    $result_to_fix->input_a_id($input_a_id);
    $result_to_fix->input_b_id($input_b_id);
    $result_to_fix->incoming_vcf_result_a($input_result_vcf_a);
    $result_to_fix->incoming_vcf_result_b($input_result_vcf_b);
    $result_to_fix->variant_type("snvs"); #to date all vcf combine operations will be snvs
    $result_to_fix->joinx_version("unknown"); # could leave this blank or something else, but I think this label works
    return 1;
}

# Get the vcf result of the non-vcf result that matches the vcf version given
sub get_vcf_result {
    my $result = shift;
    my $vcf_version = shift;
    my $aligned_reads_sample = shift;

    my $vcf_result = Genome::Model::Tools::DetectVariants2::Result::Vcf->get(
        input_id => $result->id,
        vcf_version => $vcf_version,
        aligned_reads_sample => $aligned_reads_sample,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );
    return $vcf_result;

#    print " ====================================================================\n";
#        print Data::Dumper::Dumper(@result);
#    print " --------------------------------------------------------------------\n";
#    unless(@result < 2){
#        die("Found ".scalar(@result)." vcf results for vcf_version: ".$vcf_version . " and input_id: " . $result->id);
#    }
#    my $vcf_result = (@result == 1) ? $result[0] : undef;
#    return $vcf_result;
}


1;
