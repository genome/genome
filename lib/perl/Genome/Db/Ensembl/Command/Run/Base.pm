package Genome::Db::Ensembl::Command::Run::Base;

use strict;
use warnings;
use Genome;
use Cwd;
use IO::Handle;
use File::Basename qw(dirname);

my ($VEP_DIR) = Cwd::abs_path(__FILE__) =~ /(.*)\//;
my $VEP_SCRIPT_PATH = $VEP_DIR . "/Vep.d/vep";

class Genome::Db::Ensembl::Command::Run::Base {
    is => 'Command::V2',
    doc => 'Run VEP',
    attributes_have => [
        is_vep_param => {is => 'Boolean', is_optional => 1},
    ],
    has => [
        version => {
            is => 'String',
            doc => 'version of the Variant Effects Predictor to use. Note: Version 2.7 requires using annotation version 69_37n_v3 or greater',
            valid_values => [qw(2_2 2_5 2_7)],
            is_optional => 1,
            default_value => "2_2",
        },
        input_file => {
            is => 'String',
            doc => 'File of variants to be annotated',
            is_input => 1,
            is_vep_param => 1,
        },
        format => {
            is => 'String',
            doc => 'The format of the input file, or guess to try to work out format',
            valid_values => [qw(ensembl pileup vcf hgvs id bed)],
            default_value => "bed",
            is_input => 1,
            is_vep_param => 1,
        },
        output_file => {
            is => 'String',
            doc => 'File of annotated variants.  Write to STDOUT by specifying -o STDOUT',
            is_input => 1,
            is_output => 1,
            is_vep_param => 1,
        },
        species => {
            is => 'String',
            doc => 'Species to use',
            is_optional => 1,
            default_value => 'homo_sapiens',
            is_vep_param => 1,
        },
        terms => {
            is => 'String',
            doc => 'Type of consequence terms to output',
            is_optional => 1,
            default_value => 'ensembl',
            valid_values => [qw(ensembl SO NCBI)],
            is_vep_param => 1,
        },
        sift => {
            is => 'String',
            doc => 'Add SIFT [p]rediction, [s]core or [b]oth',
            is_optional => 1,
            valid_values => [qw(p s b)],
            is_input => 1,
            is_vep_param => 1,
        },
        polyphen => {
            is => 'String',
            doc => 'Add PolyPhen [p]rediction, [s]core or [b]oth',
            is_optional => 1,
            valid_values => [qw(p s b)],
            is_input => 1,
            is_vep_param => 1,
        },
        condel => {
            is => 'String',
            doc => 'WARNING! This option is only valid for older versions of ensembl.  Use "--plugins" instead! Add Condel [p]rediction, [s]core or [b]oth.',
            is_optional => 1,
            valid_values => [qw(p s b)],
            is_input => 1,
            is_vep_param => 1,
        },
        gtf_cache => {
            is => 'Boolean',
            doc => 'Whether to use gtf_cache',
            is_optional => 1,
            default => 0,
            is_vep_param => 1,
        },
        gtf_file => {
            is => 'Text',
            doc => 'Gtf file to create cache rather than using db',
            is_optional => 1,
            is_vep_param => 1,
        },
        regulatory => {
            is => 'Boolean',
            doc => 'Look for overlap with regulatory regions.',
            default_value => 0,
            is_optional => 1,
            is_vep_param => 1,
        },
        vcf => {
            is => 'Boolean',
            doc => 'Output in VCF format instead of VEP default.',
            default_value => 0,
            is_optional => 1,
            is_vep_param => 1,
        },
        gene => {
            is => 'Boolean',
            doc => 'Force output fo Ensembl gene identifier.',
            default_value => 0,
            is_optional => 1,
            is_vep_param => 1,
        },
        most_severe => {
            is => 'Boolean',
            doc => 'Output only the most severe consequence per variation.  Transcript-specific columns will be left blank.',
            default_value => 0,
            is_optional => 1,
            is_vep_param => 1,
        },
        per_gene => {
            is => 'Boolean',
            doc => 'Output only the most severe consequence per gene.  The transcript selected is arbitrary if more than one has the same predicted consequence.',
            default_value => 0,
            is_optional => 1,
            is_vep_param => 1,
        },
        hgvs => {
            is => 'Boolean',
            doc => 'Output hgvs terms.  Must provide a reference fasta to use this option',
            default_value => 0,
            is_optional => 1,
            is_vep_param => 1,
        },
        fasta => {
            is => 'String',
            doc => 'Reference fasta file',
            is_optional => 1,
            is_vep_param => 1
        },
        hgnc => {
            is => 'Boolean',
            doc => 'NOTE: this option is ignored. The gene symbol (e.g. HGNC) will always be printed',
            is_optional => 1,
            is_input => 1,
        },
        coding_only => {
            is => 'Boolean',
            doc => 'Only return consequences that fall in the coding regions of transcripts.',
            default_value => 0,
            is_optional => 1,
            is_vep_param => 1,
        },
        force => {
            is => 'Boolean',
            doc => 'By default, the script will fail with an error if the output file already exists.  You can force the overwrite of the existing file by using this flag.',
            default_value => 0,
            is_optional => 1,
            is_vep_param => 1,
        },
        reference_build_id => {
            is => "String",
            doc => "Id of ReferenceSequence that the annotation is based on.  Only needed if --gtf-cache is specified",
            is_optional => 1,
        },
        canonical => {
            is => 'Boolean',
            doc => "Adds a flag indicating if the transcript is the canonical transcript for the gene",
            default_value => 0,
            is_optional => 1,
            is_vep_param => 1,
        },
        plugins => {
            is => 'String',
            is_many => 1,
            is_optional => 1,
            doc => 'Plugins to use.  These should be formated: PluginName@arg1@arg2,...  If you need to reference the plugin config directory in any of the args, use the placeholder PLUGIN_DIR.  For example, the Condel plugin should be called as follows: Condel@PLUGIN_DIR@b@2',
        },
        plugins_version => {
            is => 'String',
            default => "1",
            doc => 'Version of the vepplugins package to use',
        },
        quiet => {
            is => 'Boolean',
            default => 0,
            doc => "Don't print the potput of vep to the terminal",
            is_vep_param => 1,
        },
        custom => {
            is => 'String',
            is_optional => 1,
            is_many => 1,
            doc => "--custom option(s) to pass on to VEP.  Replace commas with @ symbol"
        },
    ],
    has_transient_optional => [
        _workspace => {
            is => 'Text',
            doc => "temporary to put cache and config files expected by vep",
        },
        _plugins_workspace => {
            is => 'Text',
            doc => "path to plugins config directory in workspace",
        },
    ],
};

sub help_brief {
    'Tool to run Ensembl VEP (Variant Effect Predictor)';
}

sub help_detail {
    return <<EOS
    Tool to run Ensembl VEP (Variant Effect Predictor).  For documentation on input format, see:
    http://useast.ensembl.org/info/docs/variation/vep/vep_formats.html

    It is recommended that the input file be in BED format. For example:
    5    140531    140532    T      C
    1    881906    881906    -      C
    8    12599     12602     CGT    -

    This wrapper will convert your BED to the corresponding Ensembl format, which looks like:
    5    140532    140532    T/C    +
    1    881907    881906    -/C    +
    8    12600     12602     CGT/-  +

    If using Ensembl format as input, the 5th column must always be a '+' because we always call
    variants on the forward strand. Also notice how start > stop for Ensembl format insertions.
EOS
}

sub _open_input_file {
    my $self = shift;
    my $input_file = shift;
    return Genome::Sys->open_file_for_reading($input_file) if $input_file ne '-';

    my $ioh = new IO::Handle;
    unless ($ioh->fdopen(fileno(STDIN), "r")) {
        die $self->error_message("Failed to open standard input");
    }
    return $ioh;
}

sub execute {
    my $self = shift;

    $self->resolve_format_and_input_file;
    $self->stage_plugins;
    $self->stage_cache;

    $self->run_command();

    return 1;
}

sub run_command {
    my $self = shift;
    my $cmd = $self->command;
    $self->debug_message("Running command:\n%s", $cmd);

    my %params = (
        cmd=>$cmd,
        output_files => [$self->{output_file}],
        skip_if_output_is_present => 0,
    );
    $params{input_files} = [$self->input_file] unless $self->input_file eq '-';
    $params{redirect_stdout} = '/dev/null' if $self->quiet;

    $self->api->prepend_api_path_and_execute(%params);
}

sub api {
    my $self = shift;
    return Genome::Db::Ensembl::Api->get_or_create(version => $self->ensembl_version,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME});
}

sub command {
    my $self = shift;
    return join(" ", $self->script_path, $self->string_args, $self->bool_args,
            $self->plugin_args, $self->custom_args,
            $self->db_connect_args, $self->cache_args);
}

sub workspace {
    my $self = shift;
    unless (defined $self->_workspace) {
        $self->_workspace(Genome::Sys->create_temp_directory);
    }
    return $self->_workspace;
}

sub script_path {
    my $self = shift;
    my $script_path = $self->_resolve_vep_script_path($self->ensembl_version);
    return $script_path;
}

sub ensembl_version {
    die "Must implement ensembl_version";
}

sub resolve_format_and_input_file {
    my $self = shift;
    my ($format, $input_file);
    if ($self->format eq "ensembl") {
        $input_file = $self->_verify_ensembl_input($self->input_file);
        $format = $self->format;
    } elsif ($self->format eq "bed") {
        # If bed format is input, we do a conversion to ensembl format. This is necessary
        # because ensembl has some weird conventions. (Most notably, an insertion is
        # represented by changing the start base to the end base + 1 and a deletion is represented by
        # the numbers of the nucleotides of the bases being affected:
        #
        # 1  123  122  -/ACGT
        # 1  978  980  ACT/-
        $input_file = $self->_convert_bed_to_ensembl_input($self->input_file);
        $format = "ensembl";
    } else {
        $input_file = $self->input_file;
        $format = $self->format;
    }
    $self->input_file($input_file);
    $self->format($format);
}

sub stage_plugins {
    my $self = shift;
    if ($self->plugins) {
        for my $plugin($self->plugins) {
            my $p = Genome::Db::Ensembl::VepPlugin->create(
                descriptor => $plugin,
                version => $self->plugins_version,
                staging_directory => $self->plugins_workspace,
            );
            $p->stage;
        }
    }
}

sub stage_cache {
    my $self = shift;
    $self->cache->stage($self->workspace);
}

sub plugins_workspace {
    my $self = shift;
    unless (defined $self->_plugins_workspace) {
        if ($self->plugins) {
            my $temp_plugins_dir = File::Spec->join($self->workspace, "Plugins");
            Genome::Sys->create_directory($temp_plugins_dir);
            $self->_plugins_workspace($temp_plugins_dir);
        }
    }
    return $self->_plugins_workspace;
}

sub string_args {
    my $self = shift;
    my $meta = $self->__meta__;
    
    my @all_string_args = $meta->properties(
        is_vep_param => 1,
        data_type => 'String');

    my $string_args = join( ' ',
        map {
            my $name = $_->property_name;
            my $value = $self->$name;
            defined($value) ? ("--".($name)." ".$value) : ()
        } @all_string_args
    );
    # the vep script does not understand --input_file - as read from stdin.
    # instead, we must leave out the --input_file arg to get that behavior.
    my $input_file_arg = $self->input_file eq '-' ? "" : sprintf("--input_file %s", $self->input_file);
    $string_args =~ s/--input_file ([^\s]+)/$input_file_arg/;

    return $string_args;
}

sub bool_args {
    my $self = shift;
    my $meta = $self->__meta__;
    my @all_bool_args = $meta->properties(
        is_vep_param => 1,
        data_type => 'Boolean');

    my @extra_args;
    
    if ($self->ensembl_version >= 74) {
        push @extra_args, "symbol";
    }
    else {
        push @extra_args, "hgnc";
    }

    my $bool_args = "";
    my @args_strings = map {
        my $name = $_->property_name;
        my $value = $self->$name;
        $value ? ("--".($name)) : ()
    } @all_bool_args;

    my @extra_args_strings = map {
        "--".$_;
    } @extra_args;

    $bool_args = join (' ',
        @args_strings, @extra_args_strings
    );

    return $bool_args;
}

sub custom_args {
    my $self = shift;
    my @custom_args;
    for my $custom ($self->custom) {
        my @parts = split "@", $custom;
        push @custom_args, "--custom ".join(",", @parts);
    }
    return join(" ", @custom_args);
}

sub plugin_args {
    my $self = shift;
    my @plugin_strings;
    if ($self->plugins) {
        for my $plugin($self->plugins) {
            my $p = Genome::Db::Ensembl::VepPlugin->create(
                descriptor => $plugin,
                version => $self->plugins_version,
                staging_directory => $self->plugins_workspace,
            );
            push @plugin_strings, $p->command_line_args;
        }
    }
    return join(" ", @plugin_strings);
}

sub cache_args {
    my $self = shift;
    my $cache_args = "";

    my $cache_result = $self->cache;
    if ($cache_result) {
        $self->debug_message("Using VEP cache result ".$cache_result->id);
        $cache_args = "--cache --offline --dir ".$self->workspace."/";
    }
    else {
        $self->warning_message("No cache result available, running from database");
    }
    return $cache_args;
}

sub db_connect_args {
    my $host_param = defined $ENV{GENOME_DB_ENSEMBL_HOST} ? "--host ".$ENV{GENOME_DB_ENSEMBL_HOST} : "";
    my $user_param = defined $ENV{GENOME_DB_ENSEMBL_USER} ? "--user ".$ENV{GENOME_DB_ENSEMBL_USER} : "";
    my $password_param = defined $ENV{GENOME_DB_ENSEMBL_PASS} ? "--password ".$ENV{GENOME_DB_ENSEMBL_PASS} : "";
    my $port_param = defined $ENV{GENOME_DB_ENSEMBL_PORT} ? "--port ".$ENV{GENOME_DB_ENSEMBL_PORT} : "";
    return join(" ", $host_param, $user_param, $password_param, $port_param);
}

sub _species_lookup {
    my ($self, $species) = @_;

    if ($species eq "homo_sapiens") {
        return "human";
    }
    else {
        return $species;
    }
}

sub _resolve_vep_script_path {
    my $self = shift;
    my $version = shift;

    my $api = Genome::Db::Ensembl::Api->get_or_create(version => $version);
    if (defined $api) {
        my $script_path = $api->vep_script("variant_effect_predictor.pl");
        if (-s $script_path) {
            return $script_path;
        }
        else {
            $self->warning_message("Ensembl api ".$api->id." did not have VEP script");
        }
    }
    else {
        $self->warning_message("Could not find ensembl api version ".$version);
    }
    return $VEP_SCRIPT_PATH.$self->{version}.".pl";
}

sub _verify_ensembl_input {
    my $self = shift;
    my $input_file = shift;
    my $inFh = $self->_open_input_file($input_file);
    my $tfh;
    my $tmpfile;

    # If we are reading from stdin, we won't be able to reopen the input, so dump
    # to a temp file while verifying the format
    if ($input_file eq '-') {
        ($tfh, $tmpfile) = Genome::Sys->create_temp_file;
    }

    while( my $line = $inFh->getline )
    {
        $tfh->print($line) if $tfh;

        chomp($line);
        my @F = split("\t",$line);

        #skip headers and blank lines
        next if $line =~/^#/;
        next if $line =~/^Chr/;
        next if $line =~/^$/;

        my @vars = split("/",$F[3]);
        #check SNVs
        if(($vars[0] =~ /^\w$/) && ($vars[1] =~ /^\w$/)){
            unless ($F[1] == $F[2]){
                die ("Ensembl variant format is 1-based. This line doesn't appear valid:\n$line\n");
            }
        }
        #indel insertion
        elsif(($vars[0] =~ /^-$/) && ($vars[1] =~ /\w+/)){
            unless ($F[1] == $F[2]+1){
                die ("This insertion is not in valid Ensembl format:\n$line\n");
            }
        }
        #indel deletion
        elsif(($vars[0] =~ /\w+/) && ($vars[1] =~ /^-$/)){
            unless ($F[1]+length($vars[0])-1 == $F[2]){
                die ("This deletion is not in valid Ensembl format:\n$line\n");
            }
        }
        else{
            die ("This variant is not in valid Ensembl format:\n$line\n");
        }
    }
    close($inFh);

    $input_file = $tmpfile if $tmpfile;

    return $input_file;
}

sub _convert_bed_to_ensembl_input {
    my $self = shift;
    my $input_file = shift;

    #create a tmp file for ensembl file
    my ($tfh,$tmpfile) = Genome::Sys->create_temp_file;
    unless($tfh) {
        die $self->error_message("Unable to create temporary file $!");
    }

    #convert the bed file
    my $inFh = $self->_open_input_file($self->input_file);
    while( my $line = $inFh->getline ){
        chomp($line);
        my @F = split("\t",$line);

        #skip headers and blank lines
        next if $line =~/^#/;
        next if $line =~/^Chr/;
        next if $line =~/^$/;

        #accept ref/var alleles as either slash or tab sep (A/C or A\tC)
        my @vars;
        my @suffix;
        if($F[3] =~ /\//){
            @vars = split(/\//,$F[3]);
            @suffix = @F[4..(@F-1)]
        }
        else {
            @vars = @F[3..4];
            @suffix = @F[5..(@F-1)]
        }
        $vars[0] =~ s/\*/-/g;
        $vars[0] =~ s/0/-/g;
        $vars[1] =~ s/\*/-/g;
        $vars[1] =~ s/0/-/g;

        #check SNVs
        if(($vars[0] =~ /^\w$/) && ($vars[1] =~ /^\w$/)){
            unless ($F[1] == $F[2]-1){
                die ("BED variant format is 0-based. This line doesn't appear valid:\n$line\n");
            }
            $F[1]++;
        }
        #indel insertion
        elsif(($vars[0] =~ /^-$/) && ($vars[1] =~ /\w+/)){
            unless ($F[1] == $F[2]){
                die ("This insertion is not in valid BED format:\n$line\n");
            }
            #increment the start position
            $F[1]++;
        }
        #indel deletion
        elsif(($vars[0] =~ /\w+/) && ($vars[1] =~ /^-$/)){
            unless ($F[1]+length($vars[0]) == $F[2]){
                die ("This deletion is not in valid BED format:\n$line\n");
            }
            #increment the start position
            $F[1]++;
        }
        else {
            die ("This variant is not in valid BED format:\n$line\n");
        }
        $tfh->print(join("\t",(@F[0..2],join("/",@vars),"+",@suffix)) . "\n");
    }
    close($inFh);
    return $tmpfile;
}

sub cache {
    my $self = shift;

    my $cache_result;
    my %cache_result_params;
    $cache_result_params{version} = $self->ensembl_version;
    $cache_result_params{species} = $self->_species_lookup($self->species);
    $cache_result_params{test_name} = $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME};
    if ($self->gtf_cache) {
        $cache_result_params{reference_build_id} = $self->reference_build_id;
        if (defined $self->gtf_file) {
            $cache_result_params{gtf_file_path} = $self->gtf_file;
            $cache_result_params{vep_version} = $self->version;
            $cache_result = Genome::Db::Ensembl::GtfCache->get_or_create(%cache_result_params);
        }
        else {
            $cache_result = Genome::Db::Ensembl::GtfCache->get(%cache_result_params);
        }
    }
    else {
        if ($self->sift or $self->polyphen) {
            $cache_result_params{sift} = 1;
        }
        else {
            $cache_result_params{sift} = 0;
        }
        eval {$cache_result = Genome::Db::Ensembl::VepCache->get_or_create(%cache_result_params);
        };
    }

    return $cache_result;
}

1;
