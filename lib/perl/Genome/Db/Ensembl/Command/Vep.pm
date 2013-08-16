package Genome::Db::Ensembl::Command::Vep;

use strict;
use warnings;
use Genome;
use Cwd;
use IO::Handle;
use File::Basename;

my ($VEP_DIR) = Cwd::abs_path(__FILE__) =~ /(.*)\//;
my $VEP_SCRIPT_PATH = $VEP_DIR . "/Vep.d/vep";

class Genome::Db::Ensembl::Command::Vep {
    is => 'Command::V2',
    doc => 'Run VEP',
    has => [
        version => {
            is => 'String',
            doc => 'version of the Variant Effects Predictor to use',
            valid_values => [qw(2_2 2_5 2_7)],
            is_optional => 1,
            default_value => "2_2",
        },
        input_file => {
            is => 'String',
            doc => 'File of variants to be annotated',
            is_input => 1,
        },
        format => {
            is => 'String',
            doc => 'The format of the input file, or guess to try to work out format',
            valid_values => [qw(ensembl pileup vcf hgvs id bed)],
            default_value => "bed",
            is_input => 1,
        },
        output_file => {
            is => 'String',
            doc => 'File of annotated variants.  Write to STDOUT by specifying -o STDOUT',
            is_input => 1,
            is_output => 1,
        },
        species => {
            is => 'String',
            doc => 'Species to use',
            is_optional => 1,
            default_value => 'homo_sapiens',
        },
        terms => {
            is => 'String',
            doc => 'Type of consequence terms to output',
            is_optional => 1,
            default_value => 'ensembl',
            valid_values => [qw(ensembl SO NCBI)],
        },
        sift => {
            is => 'String',
            doc => 'Add SIFT [p]rediction, [s]core or [b]oth',
            is_optional => 1,
            valid_values => [qw(p s b)],
            is_input => 1,
        },
        polyphen => {
            is => 'String',
            doc => 'Add PolyPhen [p]rediction, [s]core or [b]oth',
            is_optional => 1,
            valid_values => [qw(p s b)],
            is_input => 1,
        },
        condel => {
            is => 'String',
            doc => 'Add Condel [p]rediction, [s]core or [b]oth',
            is_optional => 1,
            valid_values => [qw(p s b)],
            is_input => 1,
        },
        gtf_cache => {
            is => 'Boolean',
            doc => 'Whether to use gtf_cache',
            is_optional => 1,
            default => 0,
        },
        gtf_file => {
            is => 'Text',
            doc => 'Gtf file to create cache rather than using db',
            is_optional => 1,
        },
        regulatory => {
            is => 'Boolean',
            doc => 'Look for overlap with regulatory regions.',
            default_value => 0,
            is_optional => 1,
        },
        vcf => {
            is => 'Boolean',
            doc => 'Output in VCF format instead of VEP default.',
            default_value => 0,
            is_optional => 1,
        },
        gene => {
            is => 'Boolean',
            doc => 'Force output fo Ensembl gene identifier.',
            default_value => 0,
            is_optional => 1,
        },
        most_severe => {
            is => 'Boolean',
            doc => 'Output only the most severe consequence per variation.  Transcript-specific columns will be left blank.',
            default_value => 0,
            is_optional => 1,
        },
        per_gene => {
            is => 'Boolean',
            doc => 'Output only the most severe consequence per gene.  The transcript selected is arbitrary if more than one has the same predicted consequence.',
            default_value => 0,
            is_optional => 1,
        },
        hgnc => {
            is => 'Boolean',
            doc => 'Adds the HGNC gene identifier (where available) to the output.',
            default_value => 0,
            is_optional => 1,
            is_input => 1,
        },
        coding_only => {
            is => 'Boolean',
            doc => 'Only return consequences that fall in the coding regions of transcripts.',
            default_value => 0,
            is_optional => 1,
        },
        force => {
            is => 'Boolean',
            doc => 'By default, the script will fail with an error if the output file already exists.  You can force the overwrite of the existing file by using this flag.',
            default_value => 0,
            is_optional => 1,
        },
        ensembl_annotation_build_id => {
            is => 'String',
            doc => "ID of ImportedAnnotation build with the desired ensembl version. \n  Current default build is: $ENV{GENOME_DB_ENSEMBL_DEFAULT_IMPORTED_ANNOTATION_BUILD}, \n  Ensembl 67_37l_v2 is: 124434505)",
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
        },
        plugins => {
            is => 'String',
            is_many => 1,
            is_optional => 1,
            doc => 'Plugins to use.  These should be formated: PluginName,arg1,arg2,...  If you need to reference the plugin config directory in any of the args, use the placeholder PLUGIN_DIR.  For example, the Condel plugin should be called as follows: Condel,PLUGIN_DIR,b,2',
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
    my $input_file = $self->input_file;
    return Genome::Sys->open_file_for_reading($input_file) if $input_file ne '-';

    my $ioh = new IO::Handle;
    unless ($ioh->fdopen(fileno(STDIN), "r")) {
        die $self->error_message("Failed to open standard input");
    }
    return $ioh;
}

sub execute {
    my $self = shift;
    # check for imported annotation build
    unless($self->ensembl_annotation_build_id) {
        $self->error_message("No ensembl annotation build specified");
        return;
    }
    my $annotation_build = Genome::Model::Build::ImportedAnnotation->get($self->ensembl_annotation_build_id);

    unless ($annotation_build) {
        $self->error_message("Could not find ImportedAnnotation build with id ".$self->ensembl_annotation_build_id);
        return;
    }

    unless ($annotation_build->get_api_paths) {
        $self->error_message("Could not find ensembl api in ImportedAnnotation build with id ".$annotation_build->id);
    }

    my $format = $self->format;
    my $input_file= $self->input_file;

    # sanity check for ensembl input, since some of us are unable to
    # keep our zero and one-based annotation files straight, and VEP
    # fails cryptically when given one-based indels

    if ($format eq "ensembl"){
        my $inFh = $self->_open_input_file;
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
    }

    # If bed format is input, we do a conversion to ensembl format. This is necessary
    # because ensembl has some weird conventions. (Most notably, an insertion is
    # represented by changing the start base to the end base + 1 and a deletion is represented by
    # the numbers of the nucleotides of the bases being affected:
    #
    # 1  123  122  -/ACGT
    # 1  978  980  ACT/-

    if ($format eq "bed"){
        #create a tmp file for ensembl file
        my ($tfh,$tmpfile) = Genome::Sys->create_temp_file;
        unless($tfh) {
            die $self->error_message("Unable to create temporary file $!");
        }

        #convert the bed file
        my $inFh = $self->_open_input_file;
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

        $format = "ensembl";
        $input_file = $tmpfile;
    }

    my $script_path = $VEP_SCRIPT_PATH.$self->{version}.".pl";
    my $string_args = "";

    #UR magic to get the string and boolean property lists
    my $meta = $self->__meta__;
    my @all_bool_args = $meta->properties(
        class_name => __PACKAGE__,
        data_type => 'Boolean');
    my @all_string_args = $meta->properties(
        class_name => __PACKAGE__,
        data_type => 'String');

    my $count = 0;
    foreach my $arg (@all_string_args) {
        if ($arg->property_name eq 'version') {
            splice @all_string_args, $count, 1;
            last;
        }
        $count++;
    }

    $count = 0;
    foreach my $arg (@all_string_args) {
        if ($arg->property_name eq 'ensembl_annotation_build_id') {
            splice @all_string_args, $count, 1;
            last;
        }
        $count++;
    }

    $count = 0;
    foreach my $arg (@all_string_args) {
        if ($arg->property_name eq 'plugins') {
            splice @all_string_args, $count, 1;
            last;
        }
        $count++;
    }

    $count = 0;
    foreach my $arg (@all_string_args) {
        if ($arg->property_name eq 'plugins_version') {
            splice @all_string_args, $count, 1;
            last;
        }
        $count++;
    }

    $count = 0;
    foreach my $arg (@all_string_args) {
        if ($arg->property_name eq 'gtf_file') {
            splice @all_string_args, $count, 1;
            last;
        }
        $count++;
    }

    $count = 0;
    foreach my $arg (@all_string_args) {
        if ($arg->property_name eq 'reference_build_id') {
            splice @all_string_args, $count, 1;
            last;
        }
        $count++;
    }


    $string_args = join( ' ',
        map {
            my $name = $_->property_name;
            my $value = $self->$name;
            defined($value) ? ("--".($name)." ".$value) : ()
        } @all_string_args
    );

    #have to replace these arg, because it may have changed (from bed -> ensembl)
    $string_args =~ s/--format (\w+)/--format $format/;
    # the vep script does not understand --input_file - as read from stdin.
    # instead, we must leave out the --input_file arg to get that behavior.
    my $input_file_arg = $input_file eq '-' ? "" : "--input_file $input_file";
    $string_args =~ s/--input_file ([^\s]+)/$input_file_arg/;

    $count = 0;
    foreach my $arg (@all_bool_args) {
        if ($arg->property_name eq 'gtf_cache') {
            splice @all_bool_args, $count, 1;
            last;
        }
        $count++;
    }


    my $bool_args = "";
    $bool_args = join (' ',
        map {
            my $name = $_->property_name;
            my $value = $self->$name;
            $value ? ("--".($name)) : ()
        } @all_bool_args
    );

    my $temp_config_dir = Genome::Sys->create_temp_directory;
    my $plugin_args = "";
    if ($self->plugins) {
        my $temp_plugins_dir = "$temp_config_dir/Plugins";
        my $temp_plugins_config_dir = "$temp_plugins_dir/config";
        my $plugins_version = $self->plugins_version;
        my $plugins_source_dir = "/usr/lib/vepplugins-$plugins_version/";
        Genome::Sys->create_directory($temp_plugins_dir);
        Genome::Sys->create_directory($temp_plugins_config_dir);
        foreach my $plugin ($self->plugins) {
            my @plugin_fields = split /\@/, $plugin;
            my $plugin_name = $plugin_fields[0];
            my $plugin_source_file = "$plugins_source_dir/$plugin_name.pm";
            if (-e $plugin_source_file){
                Genome::Sys->copy_file($plugin_source_file, "$temp_plugins_dir/$plugin_name.pm");
            }
            my $plugin_dir = "$plugins_source_dir/config/$plugin_name";
            my $temp_plugin_config_dir = "$temp_plugins_config_dir/$plugin_name/config";
            if (-d $plugin_dir) {
                `cp -r $plugin_dir $temp_plugins_config_dir`;
                foreach my $config_file (glob "$temp_plugin_config_dir/*") {
                    my $sed_cmd = "s|path/to/config/|$temp_plugins_config_dir/|";
                    `sed "$sed_cmd" $config_file > $config_file.new; mv $config_file.new $config_file`;
                }
            }
            $plugin =~ s|PLUGIN_DIR|$temp_plugin_config_dir|;
            $plugin_args .= " --plugin $plugin";
            $plugin_args =~ s/\@/,/g;
        }
    }

    my $host_param = defined $ENV{GENOME_DB_ENSEMBL_HOST} ? "--host ".$ENV{GENOME_DB_ENSEMBL_HOST} : "";
    my $user_param = defined $ENV{GENOME_DB_ENSEMBL_USER} ? "--user ".$ENV{GENOME_DB_ENSEMBL_USER} : "";
    my $password_param = defined $ENV{GENOME_DB_ENSEMBL_PASS} ? "--password ".$ENV{GENOME_DB_ENSEMBL_PASS} : "";
    my $port_param = defined $ENV{GENOME_DB_ENSEMBL_PORT} ? "--port ".$ENV{GENOME_DB_ENSEMBL_PORT} : "";
    my $cache_param;
    my $ensembl_version = Genome::Db::Ensembl::Command::Import::Run->ensembl_version_string($annotation_build->version);

    my $cache_result;
    my %cache_result_params;
    $cache_result_params{version} = $ensembl_version;
    my $rename_species = $self->_species_lookup($self->species);
    if ($rename_species) {
        $cache_result_params{species} = $rename_species;
    }
    else {
        $cache_result_params{species} = $self->species;
    }
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

    my $cmd = "$script_path $string_args $bool_args $plugin_args $host_param $user_param $password_param $port_param";

    if ($cache_result) {
        $self->status_message("Using VEP cache result ".$cache_result->id);
        $cmd = "$cmd --cache --dir ".$temp_config_dir."/";
        foreach my $file (glob $cache_result->output_dir."/*"){
            `ln -s $file $temp_config_dir`;
        }
        if ($self->gtf_cache) {
            $cmd = "$cmd --offline";
        }
    }
    else {
        $self->warning_message("No cache result available, running from database");
    }

    $self->status_message("Running command:\n$cmd");

    my %params = (
        cmd=>$cmd,
        output_files => [$self->{output_file}],
        skip_if_output_is_present => 0,
    );
    $params{input_files} = [$input_file] unless $input_file eq '-';
    $params{redirect_stdout} = '/dev/null' if $self->quiet;

    $annotation_build->prepend_api_path_and_execute(
        %params
    );
    return 1;
}

sub _species_lookup {
    my ($self, $species) = @_;

    if ($species eq "homo_sapiens") {
        return "human";
    }
    else {
        return;
    }
}


1;
