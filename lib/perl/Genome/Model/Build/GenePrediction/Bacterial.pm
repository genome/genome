package Genome::Model::Build::GenePrediction::Bacterial;

use strict;
use warnings;

use Genome;
use YAML;
use Carp;
use File::Basename;
use File::Find 'find';

class Genome::Model::Build::GenePrediction::Bacterial {
    is => 'Genome::Model::Build::GenePrediction',
    has_optional => [
        locus_suffix => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'locus_suffix' ],
        },
        run_type => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'run_type' ],
        },
    ],
};

sub locus_tag {
    my $self = shift;
    my $model = $self->model;
    return $self->locus_id . $self->run_type;
}

sub locus_id {
    my $self = shift;
    my $model = $self->model;

    my $locus_id = $model->locus_id;
    $locus_id .= $self->locus_suffix if $self->locus_suffix;
    return $locus_id;
}

sub post_allocation_initialization {
    my $self = shift;

    my($auto_suffix) = $self->input(name => 'auto_suffix');
    if($auto_suffix and $auto_suffix->value_id eq 1 and not $self->locus_suffix) {
        $self->add_input(name => 'locus_suffix', value_id => time(), value_class_name => 'UR::Value::Text'); #pick a unique suffix for this self
    }

    return 1;
}

sub assembly_name {
    my $self = shift;
    my $model = $self->model;
    return ucfirst($model->organism_name) . '_' . $self->locus_tag . '.newb.amgap';
}

sub org_dirname {
    my $self = shift;
    my $model = $self->model;
    return substr(ucfirst($model->organism_name), 0, 1) .
           substr($model->organism_name, index($model->organism_name, "_"));
}

sub organism_name {
    my $self = shift;
    my $model = $self->model;
    return ucfirst($model->organism_name);
}

sub config_file_path {
    my $self = shift;
    return $self->data_directory . '/' . $self->config_file_name;
}

sub config_file_name {
    return "config.yaml";
}

sub sequence_file_directory {
    my $self = shift;
    my $model = $self->model;
    my ($name, $path) = fileparse($model->assembly_contigs_file);
    return $path;
}

sub sequence_file_name {
    my $self = shift;
    my $model = $self->model;
    my ($name, $path) = fileparse($model->assembly_contigs_file);
    return $name;
}

sub create_config_file {
    my $self = shift;
    my $model = $self->model;
    my $config_file_path = $self->config_file_path;

    if (-e $config_file_path) {
        $self->debug_message("Removing existing configuration file at $config_file_path");
        my $unlink_rv = unlink $config_file_path;
        confess "Trouble removing configuration file at $config_file_path!" unless $unlink_rv;
    }

    $self->debug_message("Creating configuration file at $config_file_path");


    my $cell_type = uc($model->domain);
    $cell_type =~ s/(ARCHAEA|BACTERIA)L/$1/;

    my %params = (
        acedb_version    => $model->acedb_version,
        assembly_name    => $self->assembly_name,
        assembly_version => $model->assembly_version,
        cell_type        => $cell_type,
        gram_stain       => $model->gram_stain,
        locus_id         => $self->locus_id,
        locus_tag        => $self->locus_tag,
        minimum_length   => $model->minimum_sequence_length,
        ncbi_taxonomy_id => $model->ncbi_taxonomy_id,
        nr_db            => $model->nr_database_location,
        org_dirname      => $self->org_dirname,
        organism_name    => $self->organism_name,
        path             => $self->data_directory,
        pipe_version     => $model->pipeline_version,
        project_type     => $model->project_type,
        runner_count     => $model->runner_count,
        seq_file_dir     => $self->sequence_file_directory,
        seq_file_name    => $self->sequence_file_name,
        skip_acedb_parse => $model->skip_acedb_parse,
        workflowxml      => __FILE__.'.noblastp.outer.xml',
        ber_base_directory => '/gscmnt/gc9002/info/annotation/BER/autoannotate_v2.5' #FIXME This should be a parameter!
    );

    my $rv = YAML::DumpFile($config_file_path, \%params);
    unless ($rv) {
        $self->error_message("Could not create config file at $config_file_path!");
        return;
    }

    return $config_file_path;
}

sub files_ignored_by_diff {
    return qw(
        debug_file_\d+.txt
        job_err_\d+.txt
        job_\d+.txt
        debug.err.\d+

        dump_file.out
        psortb_raw_output_.+?.out
        merged.raw.bz2

        merge_valid
        prediction_valid

        rrna_screen_bsub.+(?:err|out)

        README
    );
}

sub regex_files_for_diff {
    return qw(
        phase_(\d)_ssid_\d+.ace$
        sqlite-[^-]+-\d{4}_\d{2}_\d{2}.dat$
    );
}

sub regex_for_custom_diff {
    return (
        ace         => '\.ace$',
        ks_bz2      => '\.ks\.bz2$',
        sorted_bz2  => '\.sorted\.bz2$',
        fasta       => '\.fasta$',
        sqlite      => 'sqlite-[^-]+-\d{4}_\d{2}_\d{2}.dat$',
    ), $_[0]->SUPER::regex_for_custom_diff(@_);
}

sub diff_ace {
    my ($self, $first_file, $second_file) = @_;

    my ($first_tag) = $first_file =~ m{([^/_]+)\_[^/]+$};
    my ($second_tag) = $second_file =~ m{([^/_]+)\_[^/]+$};
    my $first_md5 = `grep -v '$first_tag' $first_file | grep -v 'Database' | grep -v 'Data Version' | sort | md5sum`;
    my $second_md5 = `grep -v '$second_tag' $second_file | grep -v 'Database' | grep -v 'Data Version' | sort | md5sum`;
    return 1 if $first_md5 eq $second_md5;
    return 0;
}

sub diff_ks_bz2 {
    my ($self, $first_file, $second_file) = @_;

    my $first_md5 = `bunzip2 -dc $first_file | cut -f 1-2,4- | md5sum`;
    my $second_md5 = `bunzip2 -dc $second_file | cut -f 1-2,4- | md5sum`;
    return 1 if $first_md5 eq $second_md5;
    return 0;
}

sub diff_sorted_bz2 {
    my ($self, $first_file, $second_file) = @_;

    my $first_md5 = `bunzip2 -dc $first_file | cut -f 2- | grep -vP '\\d{2}-\\w{3}-\\d{4}' | md5sum`;
    my $second_md5 = `bunzip2 -dc $second_file | cut -f 2- | grep -vP '\\d{2}-\\w{3}-\\d{4}' | md5sum`;
    return 1 if $first_md5 eq $second_md5;
    return 0;
}

sub diff_fasta {
    my ($self, $first_file, $second_file) = @_;

    #The contig names contain the locus_tag, so cut them out
    my $first_md5 = qx(grep -v '>' $first_file | md5sum);
    my $second_md5 = qx(grep -v '>' $second_file | md5sum);
    return ($first_md5 eq $second_md5? 1 : 0);
}

sub diff_sqlite {
    my ($self, $first_file, $second_file) = @_;

    my $first_md5 = qx(cut -f 2- $first_file | md5sum);
    my $second_md5 = qx(cut -f 2- $second_file | md5sum);
    return ($first_md5 eq $second_md5? 1 : 0);
}

sub files_in_data_directory {
    my $self = shift;

    my $orgdir = join('/', $self->data_directory, $self->org_dirname);
    my $locusdir_pattern = $orgdir . '/*' . $self->locus_tag . '*';
    my ($locusdir) = glob($locusdir_pattern);
    unless ($locusdir) {
        # This directory is intermittently showing up; NFS cache problem?
        $self->error_message("Could not find locus directory: $locusdir_pattern");
        system(qq(ls "$orgdir"));
        die;
    }

    my @files;
    find({
        wanted => sub {
            my $file = $File::Find::name;
            push @files, $file;
        },
        follow => 1,
        follow_skip => 2, },
        $locusdir
    );
    return \@files;
}

sub compare_output {
    my ($self, $other_build_id) = @_;
    my $build_id = $self->build_id;
    confess "Require build ID argument!" unless defined $other_build_id;
    my $other_build = Genome::Model::Build->get($other_build_id);
    confess "Could not get build $other_build_id!" unless $other_build;

    unless ($self->model_id eq $other_build->model_id) {
        confess "Builds $build_id and $other_build_id are not from the same model!";
    }
    unless ($self->class eq $other_build->class) {
        confess "Builds $build_id and $other_build_id are not the same type!";
    }

    # Create hashes for each build, keys are paths relative to build directory and
    # values are full file paths
    my (%file_paths, %other_file_paths);
    require Cwd;
    for my $file (@{$self->files_in_data_directory}) {
        my $abs_path = Cwd::abs_path($file);
        next unless $abs_path; # abs_path returns undef if a subdirectory of file does not exist
        $file_paths{$self->full_path_to_relative($file)} = $abs_path;
    }
    for my $other_file (@{$other_build->files_in_data_directory}) {
        $other_file_paths{$other_build->full_path_to_relative($other_file)} = Cwd::abs_path($other_file);
    }

    # Now cycle through files in this build's data directory and compare with
    # corresponding files in other build's dir
    my %diffs;
    FILE: for my $rel_path (sort keys %file_paths) {
        my $abs_path = delete $file_paths{$rel_path};

        next if grep { $rel_path =~ /$_/ } $self->files_ignored_by_diff;

        warn "first abs_path ($abs_path) does not exist\n" unless (-e $abs_path);
        my $dir = $self->full_path_to_relative(dirname($abs_path));

        next FILE if -d $abs_path;

        my $regex;
        REGEX: for my $candidate ($self->regex_files_for_diff) {
            if($rel_path =~ /$candidate/) {
                $regex = $candidate;

                #check for captures to narrow the search for the matching file
                if($regex =~ /\([^?].*?\)/) {
                    my $modified_regex = $regex;
                    for($rel_path =~ /$regex/) {
                        #replace captures with their found values from the original file
                        $modified_regex =~ s/\([^?].*?\)/$_/;
                    }

                    $regex = $modified_regex;
                }

                last REGEX;
            }
        }

        # matches the same pattern
        my ($other_rel_path, $other_abs_path);
        OTHER: for my $other_file (sort keys %other_file_paths) {
            my $tag = $self->locus_tag;
            my $other_tag = $other_build->locus_tag;
            my $other_file_compare = $other_file;
            $other_file_compare =~ s/$other_tag/$tag/g;

            if($regex) {
                if($other_file_compare =~ /$regex/) {
                    $other_rel_path = $other_file;
                    $other_abs_path = delete $other_file_paths{$other_file};
                }
            } elsif($rel_path eq $other_file_compare) {
                #found
                $other_rel_path = $other_file;
                $other_abs_path = delete $other_file_paths{$other_file};
                last OTHER;
            }
        }

        # If file name doesn't match any regex, assume relative paths are the same
        unless (defined $other_rel_path and defined $other_abs_path) {
            $other_rel_path = $rel_path;
            $other_abs_path = delete $other_file_paths{$other_rel_path};
            unless (defined $other_abs_path) {
                $diffs{$rel_path} = "no file $rel_path from build $other_build_id";
                next FILE;
            }
        }

        # Check if the files end with a suffix that requires special handling. If not,
        # just do an md5sum on the files and compare
        my $diff_result = 0;
        my %matching_regex_for_custom_diff = $self->matching_regex_for_custom_diff($abs_path);
        if (keys %matching_regex_for_custom_diff > 1) {
            die "Path ($abs_path) matched multiple regex_for_custom_diff ('" . join("', '", keys %matching_regex_for_custom_diff) . "')!\n";
        }
        elsif (keys %matching_regex_for_custom_diff == 1) {
            my ($key) = keys %matching_regex_for_custom_diff;
            my $method = "diff_$key";
            unless($self->can($method)) {
                die "Custom diff method ($method) not implemented on class (" . $self->class . ").\n";
            }
            $diff_result = $self->$method($abs_path, $other_abs_path);
        }
        else {
            my $file_md5 = Genome::Sys->md5sum($abs_path);
            my $other_md5 = Genome::Sys->md5sum($other_abs_path);
            $diff_result = ($file_md5 eq $other_md5);
        }

        unless ($diff_result) {
            my $build_dir = $self->data_directory;
            my $other_build_dir = $other_build->data_directory;
            $diffs{$rel_path} = "files are not the same (diff -u {$build_dir,$other_build_dir}/$rel_path)";
        }
    }

    # Make sure the other build doesn't have any extra files
    for my $rel_path (sort keys %other_file_paths) {
        my $abs_path = delete $other_file_paths{$rel_path};

        next if grep { $rel_path =~ /$_/ } $self->files_ignored_by_diff;
        warn "second abs_path ($abs_path) does not exist\n" unless (-e $abs_path);
        my $dir = $self->full_path_to_relative(dirname($abs_path));
        next if -d $abs_path;

        $diffs{$rel_path} = "no file in build $build_id";
    }

    # Now compare metrics of both builds
    my %metric_diffs = $self->diff_metrics($other_build);
    @diffs{ keys %metric_diffs } = values %metric_diffs if %metric_diffs;

    return %diffs;
}

1;
