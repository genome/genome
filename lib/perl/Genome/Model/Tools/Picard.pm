package Genome::Model::Tools::Picard;

use strict;
use warnings;

use Genome;
use File::Basename;
use Sys::Hostname;
use Genome::Utility::AsyncFileSystem qw(on_each_line);

my $PICARD_DEFAULT = '1.46';
my $DEFAULT_MEMORY = 4;
my $DEFAULT_PERMGEN_SIZE = 64; #Mbytes
my $DEFAULT_VALIDATION_STRINGENCY = 'SILENT';

class Genome::Model::Tools::Picard {
    is  => 'Command',
    has_input => [
        use_version => { 
            is  => 'Version', 
            doc => 'Picard version to be used.',
            is_optional   => 1, 
            default_value => $PICARD_DEFAULT,
        },
        maximum_memory => {
            is => 'Integer',
            doc => 'the maximum memory (Gb) to use when running Java VM.',
            is_optional => 1,
            default_value => $DEFAULT_MEMORY,
        },
        maximum_permgen_memory => {
            is => 'Integer',
            doc => 'the maximum memory (Mbytes) to use for the "permanent generation" of the Java heap (e.g., for interned Strings)',
            is_optional => 1,
            default_value => $DEFAULT_PERMGEN_SIZE,
        },
        temp_directory => {
            is => 'String',
            doc => 'A temp directory to use when sorting or merging BAMs results in writing partial files to disk.  The default temp directory is resolved for you if not set.',
            is_optional => 1,
        },
        validation_stringency => {
            is => 'String',
            doc => 'Controls how strictly to validate a SAM file being read.',
            is_optional => 1,
            default_value => $DEFAULT_VALIDATION_STRINGENCY,
            valid_values => ['SILENT','STRICT','LENIENT'],
        },
        log_file => {
            is => 'String',
            doc => 'If provided, redirect stdout to this file.',
            is_optional => 1,
        },
        additional_jvm_options => {
            is => 'String',
            doc => 'Any additional parameters to pass to the JVM',
            is_optional => 1,
        },
        create_md5_file => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Whether to create an MD5 digest for any BAM files created.',
            default_value => 0,
        },
        #These parameters are mainly for use in pipelines
        _monitor_command => {
            is => 'Boolean',
            doc => 'Monitor output from command and warn via email if slow going',
            is_optional => 1,
            default_value => 0,
        },
        _monitor_mail_to => {
            is => 'String',
            doc => 'List of usernames to send mail to if triggered by monitor, separated by spaces (enclose list in quotes)',
            is_optional => 1,
        },
        _monitor_check_interval => {
            is => 'Integer',
            doc => 'Checks for new output each time this many seconds have elapsed',
            is_optional => 1,
        },
        _monitor_stdout_interval => {
            is => 'Integer',
            doc => 'Send the warning mail if no output detected for this many seconds',
            is_optional => 1,
        }
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run the Java toolkit Picard and work with SAM/BAM format files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt picard ...    
EOS
}

sub help_detail {                           
    return <<EOS 
More information about the Picard suite of tools can be found at http://picard.sourceforge.net/.
EOS
}

# NOTE: These are in order.
# Please put the most recent first.
my @PICARD_VERSIONS = (
    '1.77' => '/gsc/scripts/lib/java/samtools/picard-tools-1.77',
    '1.52' => '/gsc/scripts/lib/java/samtools/picard-tools-1.52',
    '1.46' => '/gsc/scripts/lib/java/samtools/picard-tools-1.46',
    '1.42' => '/gsc/scripts/lib/java/samtools/picard-tools-1.42',
    '1.40' => '/gsc/scripts/lib/java/samtools/picard-tools-1.40',
    '1.36' => '/gsc/scripts/lib/java/samtools/picard-tools-1.36',
    '1.31' => '/gsc/scripts/lib/java/samtools/picard-tools-1.31',
    '1.29' => '/gsc/scripts/lib/java/samtools/picard-tools-1.29',
    '1.25' => '/gsc/scripts/lib/java/samtools/picard-tools-1.25',
    '1.24' => '/gsc/scripts/lib/java/samtools/picard-tools-1.24',
    '1.23' => '/gsc/scripts/lib/java/samtools/picard-tools-1.23',
    'r436' => '/gsc/scripts/lib/java/samtools/picard-tools-r436', #contains a fix for when a whole library is unmapped
    '1.22' => '/gsc/scripts/lib/java/samtools/picard-tools-1.22',
    '1.21' => '/gsc/scripts/lib/java/samtools/picard-tools-1.21',
    '1.17' => '/gsc/scripts/lib/java/samtools/picard-tools-1.17',
    # old processing profiles used a different standard
    # this was supposed to be ONLY for things where we work directly from svn instead of released versions, like samtools :(
    'r116' => '/gsc/scripts/lib/java/samtools/picard-tools-1.16',
    'r107' => '/gsc/scripts/lib/java/samtools/picard-tools-1.07/',
    'r104' => '/gsc/scripts/lib/java/samtools/picard-tools-1.04/',
    'r103wu0' => '/gsc/scripts/lib/java/samtools/picard-tools-1.03/',
);

my %PICARD_VERSIONS = @PICARD_VERSIONS;

sub latest_version { $PICARD_VERSIONS[0] }

sub path_for_picard_version {
    my ($class, $version) = @_;
    $version ||= $PICARD_DEFAULT;
    my $path = $PICARD_VERSIONS{$version};
    return $path if defined $path;
    die 'No path found for picard version: '.$version;
}

sub default_picard_version {
    die "default picard version: $PICARD_DEFAULT is not valid" unless $PICARD_VERSIONS{$PICARD_DEFAULT};
    return $PICARD_DEFAULT;
}

sub picard_path {
    my $self = shift;
    return $self->path_for_picard_version($self->use_version);
}

sub run_java_vm {
    my $self = shift;
    my %params = @_;
    my $cmd = delete($params{'cmd'});
    unless ($cmd) {
        die('Must pass cmd to run_java_vm');
    }
    
    my $jvm_options = $self->additional_jvm_options || '';
    
    my $java_vm_cmd = 'java -Xmx'. $self->maximum_memory .'g -XX:MaxPermSize=' . $self->maximum_permgen_memory . 'm ' . $jvm_options . ' -cp '. $cmd;
    $java_vm_cmd .= ' VALIDATION_STRINGENCY='. $self->validation_stringency;
    $java_vm_cmd .= ' TMP_DIR='. $self->temp_directory;
    if ($self->create_md5_file) {
        $java_vm_cmd .= ' CREATE_MD5_FILE=TRUE';
    }
    if ($self->log_file) {
        $java_vm_cmd .= ' >> ' . $self->log_file;
    }
    
    $params{'cmd'} = $java_vm_cmd;
    
    if($self->_monitor_command) {
        $self->monitor_shellcmd(\%params);
    } else {
        my $result = Genome::Sys->shellcmd(%params);
    }
    
    return 1;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self->temp_directory) {
        my $base_temp_directory = Genome::Sys->base_temp_directory;
        my $temp_dir = File::Temp::tempdir($base_temp_directory .'/Picard-XXXX', CLEANUP => 1);
        Genome::Sys->create_directory($temp_dir);
        $self->temp_directory($temp_dir);
    }
    return $self;
}

sub monitor_shellcmd {
    my ($self,$shellcmd_args) = @_;
    my $check_interval = $self->_monitor_check_interval;
    my $max_stdout_interval = $self->_monitor_stdout_interval;

    unless($check_interval and $max_stdout_interval) {
        $self->error_message("Must specify a check interval ($check_interval) and a stdout interval ($max_stdout_interval) in order to monitor the command");
        die($self->error_message);
    }

    my $cmd = $shellcmd_args->{cmd};
    my $last_update = time;
    my $pid;
    my $w;
    $w = AnyEvent->timer(
        interval => $check_interval,
        cb => sub {
            if ( time - $last_update >= $max_stdout_interval) {
                my $message = <<MESSAGE;
To whom it may concern,

This command:

$cmd

Has not produced output on STDOUT in at least $max_stdout_interval seconds.

Host: %s
Perl Pid: %s 
Java Pid: %s
LSF Job: %s
User: %s

This is the last warning you will receive about this process.
MESSAGE

                undef $w;
                my $from = '"' . __PACKAGE__ . sprintf('" <%s@genome.wustl.edu>', Genome::Sys->username);

                my @to = split(' ', $self->_monitor_mail_to);
                my $to = join(', ', map { "$_\@genome.wustl.edu" } @to);
                my $subject = 'Slow ' . $self->class . ' happening right now';
                my $data = sprintf($message,
                    hostname,$$,$pid,$ENV{LSB_JOBID},Genome::Sys->username);

                my $msg = MIME::Lite->new(
                    From => $from,
                    To => $to,
                    Cc => 'apipe-run@genome.wustl.edu',
                    Subject => $subject,
                    Data => $data
                );
                $msg->send();
            }
        }
    );

    my $cv = Genome::Utility::AsyncFileSystem->shellcmd(
        %$shellcmd_args,
        '>' => on_each_line {
            $last_update = time;
            print $_[0] if defined $_[0];
        },
        '$$' => \$pid
    );
    $cv->cb(sub { undef $w });

    return $cv->recv;
}


sub parse_file_into_metrics_hashref {
    my ($class,$metrics_file,$metric_header_as_key) = @_;

    my $is_fh = Genome::Sys->open_file_for_reading($metrics_file);

    my $metric_key_index;
    unless ($metric_header_as_key) {
        if ($class->can('_metric_header_as_key')) {
            $metric_header_as_key = $class->_metric_header_as_key;
        } else {
            $class->status_message('Assuming the first column is the key for the metrics hashsref in file: '. $metrics_file);
            $metric_key_index = 0;
        }
    }
    
    my @headers;
    my %data;
    while (my $line = $is_fh->getline) {
        chomp($line);
        if ($line =~ /^## METRICS CLASS/) {
            my $next_line = $is_fh->getline;
            chomp($next_line);
            @headers = split("\t",$next_line);
            for (my $i = 0; $i < scalar(@headers); $i++) {
                unless (defined($metric_key_index)) {
                    if ($headers[$i] eq $metric_header_as_key) {
                        $metric_key_index = $i;
                    }
                }
            }
            next;
        }
        if (@headers) {
            if ($line =~ /^\s*$/) {
                last;
            } else {
                my @values = split("\t",$line);
                my $metric_key = $headers[$metric_key_index] .'-'. $values[$metric_key_index];
                for (my $i = 0; $i < scalar(@values); $i++) {
                    my $header = $headers[$i];
                    my $value = $values[$i];
                    $data{$metric_key}{$header} = $value;
                }
            }
        }
    }
    return \%data;
}

sub parse_metrics_file_into_histogram_hashref {
    my ($class,$metrics_file,$metric_header_as_key) = @_;
    
    my $as_fh = Genome::Sys->open_file_for_reading($metrics_file);
    my $metric_key_index;
    
    unless ($metric_header_as_key) {
        if ($class->can('_histogram_header_as_key')) {
            $metric_header_as_key = $class->_metric_header_as_key;
        } else {
            $class->status_message('Assuming the first column is the key for the histogram hashsref in file: '. $metrics_file);
            $metric_key_index = 0;
        }
    }

    my @headers;
    my %data;
    while (my $line = $as_fh->getline) {
        chomp($line);
        if ($line =~ /^## HISTOGRAM/) {
            my $next_line = $as_fh->getline;
            chomp($next_line);
            @headers = split("\t",$next_line);
            for (my $i = 0; $i < scalar(@headers); $i++) {
                unless (defined($metric_key_index)) {
                    if ($headers[$i] eq $metric_header_as_key) {
                        $metric_key_index = $i;
                    }
                }
            }
            next;
        }
        if (@headers) {
            if ($line =~ /^\s*$/) {
                last;
            } else {
                my $category;
                my @values = split("\t",$line);
                unless (scalar(@values) == scalar(@headers)) {
                    $DB::single=1;
                    next;
                }
                my $metric_key = $headers[$metric_key_index] .'-'. $values[$metric_key_index];
                for (my $i = 0; $i < scalar(@values); $i++) {
                    my $header = $headers[$i];
                    my $value = $values[$i];
                    $data{$metric_key}{$header} = $value;
                }
            }
        }
    }
    return \%data;
}


1;

