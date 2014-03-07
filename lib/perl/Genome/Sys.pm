package Genome::Sys;

use strict;
use warnings;
use autodie qw(chown);
use Genome;
use Cwd;
use File::Path;
use File::Spec;
use File::Basename;
use File::Copy;
use Carp;
use IO::File;
use LWP::Simple qw(getstore RC_OK);
use List::MoreUtils "each_array";
use Set::Scalar;
use Digest::MD5;
use JSON;
use Params::Validate qw(:types validate_pos);
use POSIX qw(EEXIST);

# these are optional but should load immediately when present
# until we can make the Genome::Utility::Instrumentation optional (Net::Statsd deps)
for my $opt (qw/Genome::Sys::Lock Genome::Sys::Log/) {
    eval "use $opt";
}

our $VERSION = $Genome::VERSION;

class Genome::Sys {};

sub concatenate_files {
    my ($self,$inputs,$output) = @_;
    unless (ref($inputs) and $output) {
        die 'expected \\@input_files $output_file as parameters';
    }
    my ($quoted_output_file, @quoted_input_files) = $self->quote_for_shell($output,@$inputs);
    Genome::Sys->shellcmd(cmd => "cat @quoted_input_files >$quoted_output_file");
}

sub quote_for_shell {
    # this is needed until shellcmd supports an array form,
    # which is difficult because we go to a bash sub-shell by default
    my $class = shift;
    my @quoted = @_;
    for my $value (@quoted) {
        $value =~ s|\\|\\\\|g;
        $value =~ s/\"/\\\"/g;
        $value = "\"${value}\"";
        print STDERR $value,"\n";
    }
    if (wantarray) {
        return @quoted
    }
    else {
        return $quoted[0];
    }
}


sub disk_usage_for_path {
    my $self = shift;
    my $path = shift;
    my %params = @_;

    # Normally disk_usage_for_path returns nothing if it encounters an error;
    # this allows the caller to override the default behavior, returning a
    # number even if du emits errors (for instance, because the user does not
    # have read permissions for a single subfolder out of many).
    my $allow_errors = delete $params{allow_errors};
    if (keys %params) {
        die $self->error_message("Unknown parameters for disk_usage_for_path(): " . join(', ', keys %params));
    }

    unless (-d $path) {
        $self->error_message("Path $path does not exist!");
        return;
    }

    return unless -d $path;
    my $cmd = "du -sk $path 2>&1";
    my @du_output = split( /\n/, qx{$cmd} );
    my $last_line = pop @du_output;
    my $kb_used = ( split( ' ', $last_line, 2 ) )[0];

    # if we couldn't parse the size count, print all output and return nothing
    unless (Scalar::Util::looks_like_number($kb_used)) {
        push @du_output, $last_line;
        $self->error_message("du output is not a number:\n" . join("\n", @du_output));
        return;
    }

    # if there were lines besides the size count, emit warnings and (potentially) return nothing
    if (@du_output) {
        $self->warning_message($_) for (@du_output);
        unless ($allow_errors) {
            $self->error_message("du encountered errors");
            return;
        }
    }

    return $kb_used;
}

# locking has been moved into its own module

sub lock_resource {
    shift;
    Genome::Sys::Lock->lock_resource(@_);
}

sub unlock_resource {
    shift;
    Genome::Sys::Lock->unlock_resource(@_);
}

# directory manipulation

sub validate_existing_directory {
    my ($self, $directory) = @_;

    unless ( defined $directory ) {
        Carp::croak("Can't validate_existing_directory: No directory given");
    }

    unless ( -e $directory ) {
        Carp::croak("Can't validate_existing_director: $directory: Path does not exist");
    }


    unless ( -d $directory ) {
        Carp::croak("Can't validate_existing_director: $directory: path exists but is not a directory");
    }

    return 1;
}

sub validate_directory_for_read_access {
    my ($self, $directory) = @_;

    # Both underlying methods throw their own exceptions
    $self->validate_existing_directory($directory)
        or return;

    return $self->_can_read_from_directory($directory);
}

sub validate_directory_for_write_access {
    my ($self, $directory) = @_;

    # Both underlying methods throw their own exceptions
    $self->validate_existing_directory($directory)
        or return;

    return $self->_can_write_to_directory($directory);
}

sub validate_directory_for_read_write_access {
    my ($self, $directory) = @_;

    # All three underlying methods throw their own exceptions
    $self->validate_existing_directory($directory)
        or return;

    $self->_can_read_from_directory($directory)
        or return;

    return $self->_can_write_to_directory($directory);
}

sub recursively_validate_directory_for_read_write_access {
    my ($self, $directory) = @_;

    my $wanted = sub {
        my $full_path = $File::Find::name;
        if (-f $full_path) {
            eval {
                Genome::Sys->validate_file_for_reading($full_path);
                Genome::Sys->validate_file_for_writing_overwrite($full_path);
            };

            if ($@) {
                Carp::croak "Cannot read or write $full_path in $directory!";
            }
        }
    };

    find($wanted, $directory);
    return 1;
}

sub _can_read_from_directory {
    my ($self, $directory) = @_;

    unless ( -r $directory ) {
        Carp::croak("Directory ($directory) is not readable");
    }

    return 1;
}

sub _can_write_to_directory {
    my ($self, $directory) = @_;

    unless ( -w $directory ) {
        Carp::croak("Directory ($directory) is not writable");
    }

    return 1;
}

# MD5

sub md5sum {
    my ($self, $file) = @_;

    my $digest;

    my $fh = IO::File->new($file);
    unless ($fh) {
        Carp::croak("Can't open file ($file) to md5sum: $!");
    }
    my $d = Digest::MD5->new;
    $d->addfile($fh);
    $digest = $d->hexdigest;
    $fh->close;

    return $digest;
}

sub md5sum_data {
    my ($self, $data) = @_;
    unless (defined $data) {
        Carp::croak('No data passed to md5sum_data');
    }
    my $digest = Digest::MD5->new;
    $digest->add($data);
    return $digest->hexdigest;
}
# API for accessing software and data by version

sub snapshot_revision {
    my $class = shift;

    # Previously we just used UR::Util::used_libs_perl5lib_prefix but this did not
    # "detect" a software revision when using code from PERL5LIB or compile-time
    # lib paths. Since it is common for developers to run just Genome from a Git
    # checkout we really want to record what versions of UR, Genome, and Workflow
    # were used.

    my @orig_inc = @INC;
    my @libs = ($INC{'UR.pm'}, $INC{'Genome.pm'});
    die $class->error_message('Did not find both modules loaded (UR and Genome).') unless @libs == 2;

    # assemble list of "important" libs
    @libs = map { File::Basename::dirname($_) } @libs;
    push @libs, UR::Util->used_libs;

    # remove trailing slashes
    map { $_ =~ s/\/+$// } (@libs, @orig_inc);

    @libs = $class->_uniq(@libs);

    # preserve the list order as appeared @INC
    my @inc;
    for my $inc (@orig_inc) {
        push @inc, grep { $inc eq $_ } @libs;
    }

    @inc = $class->_uniq(@inc);

    @inc = $class->_simplify_inc(@inc) if $class->can('_simplify_inc');

    return join(':', @inc);
}

sub _uniq {
    my $self = shift;
    my @list = @_;
    my %seen = ();
    my @unique = grep { ! $seen{$_} ++ } @list;
    return @unique;
}

# access to paths to code and data

sub dbpath {
    my ($class, $name, $version) = @_;
    my $envname = $class->dbname_to_envname($name);
    my $dbpath;
    if ($envname and $ENV{$envname}) {
        $dbpath = $ENV{$envname};
        print STDERR "Using '$dbpath' from $envname.\n";
    } else {
        unless ($version) {
            die "Genome::Sys dbpath must be called with a database name and a version. " .
            "Use 'latest' for the latest installed version.";
        }
        $dbpath = $class->_find_in_genome_db_paths($name, $version);
    }
    return $dbpath;
}

sub dbname_to_envname {
    my ($class, $name) = @_;
    my %envmap = (
        'genome-music-testdata' => 'GENOME_DB_MUSIC_TESTDATA',
        'cosmic' => 'GENOME_DB_COSMIC',
        'omim' => 'GENOME_DB_OMIM',
        'pfam' => 'GENOME_DB_PFAM',
    );
    return $envmap{$name};
}

sub _find_in_genome_db_paths {
    my ($class, $name, $version) = @_;

    my $base_dirs = $ENV{GENOME_DB};
    my $subdir = "$name/$version";

    my @base_dirs = split(':',$base_dirs);
    my @dirs =
        map { -l $_ ? Cwd::abs_path($_) : ($_) }
        map {
            my $path = join("/",$_,$subdir);
            (-e $path ? ($path) : ())
        }
        @base_dirs;
    return $dirs[0];
}

# renamed for consistency with a variety of sw_ methods.
*swpath = \&sw_path;

sub sw_path {
    my ($class, $pkg_name, $version, $app_name) = @_;
    $app_name ||= $pkg_name;

    unless ($version) {
        die "Genome::Sys swpath must be called with a pkg name and a version. The optional executable name defaults to the pkg_name. " .
            "Use the version 'latest' for the latest installed version.";
    }

    # check the default path for the app
    my %map = $class->sw_version_path_map($pkg_name, $app_name);
    my $path = $map{$version};
    if ($path) {
        return $path;
    }

    # older legacy software has an unversioned executable
    # this is only supported if the version passed in is "latest"
    $path = `which $app_name`;
    if ($path = `which $app_name`) {
        # unversioned install
        # see if it's a symlink to something in a versioned tree
        chomp $path;
        $path = readlink($path) while -l $path;
        if ($version eq 'latest') {
            return $path;
        }
        else {
            die $class->error_message("Failed to find $pkg_name at version $version. " .
                "The default version is at $path.");
        }
    }

    die $class->error_message("Failed to find app $app_name (package $pkg_name) at version $version!");
}

sub jar_path {
    my ($class, $jar_name, $version) = @_;
    # check the default path
    my %map = $class->jar_version_path_map($jar_name);
    my $path = $map{$version};
    if ($path) {
        return $path;
    }
    die $class->error_message("Failed to find jar $jar_name at version $version");
}

sub jar_version_path_map {
    my ($class, $pkg_name) = @_;

    my %versions;
    my @dirs = split(':', $ENV{GENOME_JAR_PATH});

    for my $dir (@dirs) {
        my $prefix = "$dir/$pkg_name-";
        my @version_paths = grep { -e $_ } glob("$prefix*");
        next unless @version_paths;
        my $prefix_len = length($prefix);
        for my $version_path (@version_paths) {
            my $version = substr($version_path,$prefix_len);
            $version =~ s/.jar$//;
            if (substr($version,0,1) eq '-') {
                $version = substr($version,1);
            }
            next unless $version =~ /[0-9\.]/;
            next if -l $version_path;
            $versions{$version} = $version_path;
        }
    }
    return %versions;
}

sub sw_version_path_map {
    my ($class, $pkg_name, $app_name) = @_;
    $app_name ||= $pkg_name;

    # find software installed as .debs with versioned packages
    # packaged software should have a versioned executable like /usr/bin/myapp1.2.3 or /usr/bin/mypackage-myapp1.2.3 in the bin.
    my %versions1;
    my @dirs1 = (split(':',$ENV{PATH}), "~/gsc-pkg-bio/");
    my @sw_ignore = split(':',$ENV{GENOME_SW_IGNORE} || '');
    for my $dir1 (@dirs1) {
        if (grep { index($dir1,$_) == 0 } @sw_ignore) {
            # skip directories starting with something in @sw_ignore
            next;
        }
        for my $prefix ("$dir1/$pkg_name-$app_name-", "$dir1/$app_name-", "$dir1/$pkg_name-$app_name", "$dir1/$app_name") {
            my @version_paths = grep { -e $_ } glob("$prefix*");
            next unless @version_paths;
            my $prefix_len = length($prefix);
            for my $version_path (@version_paths) {
                my $version = substr($version_path,$prefix_len);
                if (substr($version,0,1) eq '-') {
                    $version = substr($version,1);
                }
                next unless $version =~ /[0-9\.]/;
                if (grep { index($version_path,$_) == 0 } @sw_ignore) {
                    next;
                }
                $versions1{$version} = $version_path;
            }
        }
    }

    # find software installed under $GENOME_SW
    # The pattern for these has varied over time.  They will be like will be like:
    #   $GENOME_SW/$pkg_name/{$pkg_name-,,*}$version/{.,bin,scripts}/{$pkg_name,}{-,}$app_name{-,}{$version,}
    sub _common_prefix_length {
        my ($first,@rest) = @_;
        my $len;
        for ($len = 1; $len <= length($first); $len++) {
            for my $other (@rest) {
                if (substr($first,0,$len) ne substr($other,0,$len)) {
                    return $len-1;
                }
            }
        }
        return $len-1;
    }

    my %pkgdirs;
    my @dirs2 = split(':',$ENV{GENOME_SW});  #most of the system expects this to be one value not-colon separated currently
    for my $dir2 (@dirs2) {
        # one subdir will exist per application
        my @app_subdirs = glob("$dir2/$pkg_name");
        for my $app_subdir (@app_subdirs) {
            # one subdir under that will exist per version
            my @version_subdirs = grep { -e $_ and not -l $_ and $_ !~ /README/ and /\d/ } glob("$app_subdir/*");
            next unless @version_subdirs;

            # if some subdirectories repeat the pkg name, all of those we consider must
            my @some = grep { index(File::Basename::basename($_),$pkg_name) == 0 } @version_subdirs;
            @version_subdirs = @some if @some;

            my $len = _common_prefix_length(@version_subdirs);
            my $prefix = substr($version_subdirs[0],0,$len);
            if ($prefix =~ /(-|)(v|)([\d\.pa]+)$/) {
                $len = $len - length($3);
            }
            for my $version_subdir (@version_subdirs) {
                if (grep { index($version_subdir,$_) == 0 } @sw_ignore) {
                    next;
                }
                my $version = substr($version_subdir,$len);
                if (substr($version,0,1) =~ /[0-9]/) {
                    $pkgdirs{$version} = $version_subdir;
                }
            }
        }
    }

    # When a version number ends in -64, we are just saying it is 64-bit,
    # which is default for everything now.  We don't want to even return
    # the 32-bit only versions.
    # Trim the version number down, and have that key point to the path of the 64-bit version.
    # This will sometimes stomp on a 32-bit version directory if it exists (hopefully).
    # If 32-bit directories exist and there is no 64-bit version to stomp on them,
    # we never know for sure they are not just new installs which presume to be 64-bit.
    my @v64 = grep { /[-_]64$/ } keys %pkgdirs;
    for my $version_64 (@v64) {
        my ($version_no64) = ($version_64 =~ /^(.*?)([\.-_][^\.]+)64$/);
        $pkgdirs{$version_no64} = delete $pkgdirs{$version_64};
    }

    # now resolve the actual path to the executable $app_name
    my %versions2;
    for my $version (keys %pkgdirs) {
        my $dir = $pkgdirs{$version};
        my @subdirs = qw/. bin scripts/;
        my @basenames = ("$app_name-$version","$app_name$version",$app_name);
        if ($pkg_name ne $app_name) {
            @basenames = map { ("$pkg_name-$_", "$pkg_name$_", $_) } @basenames;
        }
        for my $basename (@basenames) {
            for my $subdir (@subdirs) {
                my $path = "$dir/$subdir/$basename";
                if (-e $path) {
                    $path =~ s|/\./|/|m;
                    $versions2{$version} = $path;
                    last;
                }
            }
            last if $versions2{$version};
        }
        unless ($versions2{$version}) {
            $class->debug_message("Found $pkg_name at version $version at $dir, but no executable (named any of: @basenames) in subdirs: @subdirs!");
        }
    }

    # prefer packaged versions over unpackaged
    my %all_versions = (%versions2,%versions1);

    # return sorted
    return map { $_ => $all_versions{$_} } sort keys %all_versions;
}

sub sw_versions {
    my ($self,$pkg_name,$app_name) = @_;
    my @map = ($self->sw_version_path_map($pkg_name,$app_name));
    my $n = 0;
    return map { $n++; ($n % 2 ? $_ : ()) } @map;
}


#####
# Temp file management
#####

sub _temp_directory_prefix {
    my $self = shift;
    my $base = join("_", map { lc($_) } split('::',$self->class));
    return $base;
}

our $base_temp_directory;
sub base_temp_directory {
    my $self = shift;
    my $class = ref($self) || $self;
    my $template = shift;

    my $id;
    if (ref($self)) {
        return $self->{base_temp_directory} if $self->{base_temp_directory};
        $id = $self->id;
    }
    else {
        # work as a class method
        return $base_temp_directory if $base_temp_directory;
        $id = '';
    }

    unless ($template) {
        my $prefix = $self->_temp_directory_prefix();
        $prefix ||= $class;
        my $time = $self->__context__->now;

        $time =~ s/[\s\: ]/_/g;
        $template = "/gm-$prefix-$time-$id-XXXX";
        $template =~ s/ /-/g;
    }

    # See if we're running under LSF and LSF gave us a directory that will be
    # auto-cleaned up when the job terminates
    my $tmp_location = $ENV{'TMPDIR'} || File::Spec->tmpdir();
    if ($ENV{'LSB_JOBID'}) {
        my $lsf_possible_tempdir = sprintf("%s/%s.tmpdir", $tmp_location, $ENV{'LSB_JOBID'});
        $tmp_location = $lsf_possible_tempdir if (-d $lsf_possible_tempdir);
    }
    # tempdir() thows its own exception if there's a problem

    # For debugging purposes, allow cleanup to be disabled
    my $cleanup = 1;
    if($ENV{'GENOME_SYS_NO_CLEANUP'}) {
        $cleanup = 0;
    }
    my $dir = File::Temp::tempdir($template, DIR=>$tmp_location, CLEANUP => $cleanup);

    $self->create_directory($dir);

    if (ref($self)) {
        return $self->{base_temp_directory} = $dir;
    }
    else {
        # work as a class method
        return $base_temp_directory = $dir;
    }

    unless ($dir) {
        Carp::croak("Unable to determine base_temp_directory");
    }

    return $dir;
}

our $anonymous_temp_file_count = 0;
sub create_temp_file_path {
    my $self = shift;
    my $name = shift;
    unless ($name) {
        $name = 'anonymous' . $anonymous_temp_file_count++;
    }
    my $dir = $self->base_temp_directory;
    my $path = $dir .'/'. $name;
    if (-e $path) {
        Carp::croak "temp path '$path' already exists!";
    }

    if (!$path or $path eq '/') {
        Carp::croak("create_temp_file_path() failed");
    }

    return $path;
}

sub create_temp_file {
    my $self = shift;
    my $path = $self->create_temp_file_path(@_);
    my $fh = IO::File->new($path, '>');
    unless ($fh) {
        Carp::croak "Failed to create temp file $path: $!";
    }
    return ($fh,$path) if wantarray;
    return $fh;
}

sub create_temp_directory {
    my $self = shift;
    my $path = $self->create_temp_file_path(@_);
    $self->create_directory($path);
    return $path;
}

#####
# Basic filesystem operations
#####

sub copy_file {
    my ($self, $file, $dest) = @_;

    $self->validate_file_for_reading($file)
        or Carp::croak("Cannot open input file ($file) for reading!");

    $self->validate_file_for_writing($dest)
        or Carp::croak("Cannot open output file ($dest) for writing!");

    # Note: since the file is validate_file_for_reading, and the dest is validate_file_for_writing,
    #  the files can never be exactly the same.

    unless ( File::Copy::copy($file, $dest) ) {
        Carp::croak("Can't copy $file to $dest: $!");
    }

    return 1;
}

sub move_file {
    my ($self, $file, $dest) = @_;

    $self->validate_file_for_reading($file)
        or Carp::croak("Cannot open input file ($file) for reading!");

    $self->validate_file_for_writing($dest)
        or Carp::croak("Cannot open output file ($dest) for writing!");

    # Note: since the file is validate_file_for_reading, and the dest is validate_file_for_writing,
    #  the files can never be exactly the same.

    unless ( File::Copy::move($file, $dest) ) {
        Carp::croak("Can't move $file to $dest: $!");
    }

    return 1;
}

sub tar {
    my ($class, %params) = @_;
    my $tar_path = delete $params{tar_path};
    my $input_directory = delete $params{input_directory};
    my $input_pattern = delete $params{input_pattern};
    $input_pattern = '*' unless defined $input_pattern;
    my $options = delete $params{options};
    $options = '-cf' unless defined $options;

    if (%params) {
        Carp::confess "Extra parameters given to tar method: " . join(', ', sort keys %params);
    }

    unless ($tar_path) {
        Carp::confess "Not given path at which tar should be created!";
    }
    if (-e $tar_path) {
        Carp::confess "File exists at $tar_path, refusing to overwrite with new tarball!";
    }

    unless ($input_directory) {
        Carp::confess "Not given directory containing input files!";
    }
    unless (-d $input_directory) {
        Carp::confess "No input directory found at $input_directory";
    }

    my $current_directory = getcwd;
    unless (chdir $input_directory) {
        Carp::confess "Could not change directory to $input_directory";
    }

    if (Genome::Sys->directory_is_empty($input_directory)) {
        Carp::confess "Cannot create tarball for empty directory $input_directory!";
    }

    my $cmd = "tar $options $tar_path $input_pattern";
    my $rv = Genome::Sys->shellcmd(
        cmd => $cmd,
    );
    unless ($rv) {
        Carp::confess "Could not create tar file at $tar_path containing files in " .
            "$input_directory matching pattern $input_pattern";
    }

    unless (chdir $current_directory) {
        Carp::confess "Could not change directory back to $current_directory";
    }
    return 1;
}

sub untar {
    my ($class, %params) = @_;
    my $tar_path = delete $params{tar_path};
    my $target_directory = delete $params{target_directory};
    my $delete_tar = delete $params{delete_tar};

    if (%params) {
        Carp::confess "Extra parameters given to untar method: " . join(', ', sort keys %params);
    }
    unless ($tar_path) {
        Carp::confess "Not given path to tar file to be untarred!";
    }
    unless (-e $tar_path) {
        Carp::confess "No file found at $tar_path!";
    }
    $target_directory = getcwd unless $target_directory;
    $delete_tar = 0 unless defined $delete_tar;

    my $current_directory = getcwd;
    unless (chdir $target_directory) {
        Carp::confess "Could not change directory to $target_directory";
    }

    my $rv = Genome::Sys->shellcmd(
        cmd => "tar -xf $tar_path",
    );
    unless ($rv) {
        Carp::confess "Could not untar $tar_path into $target_directory";
    }

    unless (chdir $current_directory) {
        Carp::confess "Could not change directory back to $current_directory";
    }

    if ($delete_tar) {
        unlink $tar_path;
    }

    return 1;
}

sub directory_is_empty {
    my ($class, $directory) = @_;
    my @files = glob("$directory/*");
    if (@files) {
        return 0;
    }
    return 1;
}

sub rsync_directory {
    my ($class, %params) = @_;
    my $source_dir = delete $params{source_directory};
    my $target_dir = delete $params{target_directory};
    my $pattern = delete $params{file_pattern};

    unless ($source_dir) {
        Carp::confess "Not given directory to copy from!";
    }
    unless (-d $source_dir) {
        Carp::confess "No directory found at $source_dir";
    }
    unless ($target_dir) {
        Carp::confess "Not given directory to copy to!";
    }
    unless (-d $target_dir) {
        Genome::Sys->create_directory($target_dir);
    }
    $pattern = '' unless $pattern;

    my $source = join('/', $source_dir, $pattern);
    my $rv = Genome::Sys->shellcmd(
        cmd => "rsync -rlHpgt $source $target_dir",
    );
    unless ($rv) {
        confess "Could not copy data matching pattern $source to $target_dir";
    }
    return 1;
}

sub line_count {
    my ($self, $path) = @_;
    my ($line_count) = qx(wc -l $path) =~ /^(\d+)/;
    return $line_count;
}

sub create_directory {
    my ($self, $directory) = @_;

    unless ( defined $directory ) {
        Carp::croak("Can't create_directory: No path given");
    }

    # FIXME do we want to throw an exception here?  What if the user expected
    # the directory to be created, not that it already existed
    return $directory if -d $directory;

    # have to set umask, make_path's mode/umask option is not sufficient
    my $umask = umask;
    umask $ENV{GENOME_SYS_UMASK};
    make_path($directory); # not from File::Path
    umask $umask;

    return $directory;
}

# File::Path::make_path says it lets you specify group but it always seemed
# to be overrided by setgid.  So we are implenting the recursive mkdir here.
sub make_path {
    my ($path) = validate_pos(@_, {type => SCALAR});

    my $gid = gidgrnam($ENV{GENOME_SYS_GROUP});

    my @dirs = File::Spec->splitdir($path);
    for (my $i = 0; $i < @dirs; $i++) {
        my $subpath = File::Spec->catdir(@dirs[0..$i]);

        my $rv = mkdir $subpath;
        my $mkdir_errno = $!;
        if ($rv) {
            chown -1, $gid, $subpath;
        } else {
            if ($mkdir_errno == EEXIST) {
                next;
            } else {
                Carp::confess("While creating path ($path), failed to create " .
                    "directory ($subpath) because ($!)");
            }
        }
    }

    unless (-d $path) {
        die "directory does not exist: $path";
    }

}

sub gidgrnam {
    my $group_name = shift;
    (getgrnam($group_name))[2];
}

sub create_symlink {
    my ($class, $target, $link) = @_;

    unless ( defined $target ) {
        Carp::croak("Can't create_symlink: no target given");
    }

    unless ( defined $link ) {
        Carp::croak("Can't create_symlink: no 'link' given");
    }

    if ( -e $link ) { # the link exists and points to spmething
        Carp::croak("Link ($link) for target ($target) already exists.");
    }

    if ( -l $link ) { # the link exists, but does not point to something
        Carp::croak("Link ($link) for target ($target) is already a link.");
    }

    unless ( symlink($target, $link) ) {
        Carp::croak("Can't create link ($link) to $target\: $!");
    }

    return 1;
}

# Return a list of strings representing the filenames/directories in a given directory.
sub list_directory {
    my ($class, $directory_name) = @_;
    opendir my($dh), $directory_name or die "Couldn't open dir '$directory_name': $!";
    my @children = readdir $dh;
    closedir $dh;

    # remove . and .. if they are in the list
    my @results = grep(!/^[\.]*$/, @children);
    return @results
}

#  create symlinks of the contents of the target_dir.
# BEFORE:
#  target_dir/foo/bar
#  target_dir/baz
#  link_dir/ (no foo or baz)
# AFTER:
#  link_dir/foo -> target_dir/foo
#  link_dir/baz -> target_dir/baz
sub symlink_directory {
    my ($class, $target_dir, $link_dir) = @_;

    if(-e $link_dir) {
        Carp::croak("The link_dir ($link_dir) exists and is not a directory!") unless(-d $link_dir);

        my $target_filenames = Set::Scalar->new(Genome::Sys->list_directory($target_dir));
        my $link_filenames = Set::Scalar->new(Genome::Sys->list_directory($link_dir));

        # check for intersection of link_filenames with target_filenames and die.
        my $intersection = $target_filenames->intersection($link_filenames);
        if(@$intersection) {
            Carp::croak("Cannot symlink directory because the following\n" .
                "are in both the target_dir and the link_dir:\n" .
                Data::Dumper::Dumper(@$intersection));
        }

        # symlink all the things in $target_dir
        my @target_fullpaths = map("$target_dir/$_", @$target_filenames);
        my @link_fullpaths = map("$link_dir/$_", @$target_filenames);
        my $ea = each_array(@target_fullpaths, @link_fullpaths);
        while( my ($target, $link) = $ea->() ) {
            Genome::Sys->create_symlink($target, $link);
        }
    } else {
            Carp::croak("The link_dir ($link_dir) doesn't exist.");
    }
}

sub create_symlink_and_log_change {
    my $class  = shift || die;
    my $owner  = shift || die;
    my $target = shift || die;
    my $link   = shift || die;

    $class->create_symlink($target, $link);

    # create a change record so that if the databse change is undone this symlink will be removed
    my $symlink_undo = sub {
        $owner->status_message("Removing symlink ($link) due to database rollback.");
        unlink $link;
    };
    my $symlink_change = UR::Context::Transaction->log_change(
        $owner, 'UR::Value', $link, 'external_change', $symlink_undo
    );
    unless ($symlink_change) {
        die $owner->error_message("Failed to log symlink change.");
    }

    return 1;
}

sub read_file {
    my ($self, $fname) = @_;
    my $fh = $self->open_file_for_reading($fname);
    Carp::croak "Failed to open file $fname! " . $self->error_message() . ": $!" unless $fh;
    if (wantarray) {
        my @lines = $fh->getlines;
        return @lines;
    }
    else {
        my $text = do { local( $/ ) ; <$fh> } ;  # slurp mode
        return $text;
    }
}

sub write_file {
    my ($self, $fname, @content) = @_;
    my $fh = $self->open_file_for_writing($fname);
    Carp::croak "Failed to open file $fname! " . $self->error_message() . ": $!" unless $fh;
    for (@content) {
        $fh->print($_) or Carp::croak "Failed to write to file $fname! $!";
    }
    $fh->close or Carp::croak "Failed to close file $fname! $!";
    return $fname;
}

sub _open_file {
    my ($self, $file, $rw) = @_;
    if ($file eq '-') {
        if ($rw eq 'r') {
            return 'STDIN';
        }
        elsif ($rw eq 'w') {
            return 'STDOUT';
        }
        else {
            die "cannot open '-' with access '$rw': r = STDIN, w = STDOUT!!!";
        }
    }
    my $fh = (defined $rw) ? IO::File->new($file, $rw) : IO::File->new($file);
    return $fh if $fh;
    Carp::croak("Can't open file ($file) with access '$rw': $!");
}

sub validate_file_for_reading {
    my ($self, $file) = @_;

    unless ( defined $file ) {
        Carp::croak("Can't validate_file_for_reading: No file given");
    }

    if ($file eq '-') {
        return 1;
    }

    unless (-e $file ) {
        Carp::croak("File ($file) does not exist");
    }

    unless (-f $file) {
        Carp::croak("File ($file) exists but is not a plain file");
    }

    unless ( -r $file ) {
        Carp::croak("Do not have READ access to file ($file)");
    }

    return 1;
}

sub validate_file_for_writing {
    my ($self, $file) = @_;

    unless ( defined $file ) {
        Carp::croak("Can't validate_file_for_writing: No file given");
    }

    if ($file eq '-') {
        return 1;
    }

    if ( -s $file ) {
        Carp::croak("Can't validate_file_for_writing: File ($file) has non-zero size, refusing to write to it");
    }

    # FIXME there is a race condition where the path could go away or become non-writable
    # between the time this method returns and the time we actually try opening the file
    # for writing

    # validate_file_for_writing_overwrite throws its own exceptions if there are problems
    return $self->validate_file_for_writing_overwrite($file);
}


sub validate_file_for_writing_overwrite {
    my ($self, $file) = @_;

    unless ( defined $file ) {
        Carp::croak("Can't validate_file_for_writing_overwrite: No file given");
    }

    my ($name, $dir) = File::Basename::fileparse($file);
    unless ( $dir ) {
        Carp::croak("Can't validate_file_for_writing_overwrite: Can't determine directory from pathname ($file)");
    }

    unless ( -w $dir ) {
        Carp::croak("Can't validate_file_for_writing_overwrite: Do not have WRITE access to directory ($dir) to create file ($name)");
    }

    # FIXME same problem with the race condition as noted at the end of validate_file_for_writing()
    return 1;
}
sub open_file_for_reading {
    my ($self, $file) = @_;

    $self->validate_file_for_reading($file)
        or return;

    # _open_file throws its own exception if it doesn't work
    return $self->_open_file($file, 'r');
}

sub download_file_to_directory {
    my ($self, $url, $destination_dir) = @_;

    unless (-d $destination_dir){
        Carp::croak("You wanted to download $url to $destination_dir but that directory doesn't exist!");
    }

    my $resp =  getstore($url, $destination_dir . "/" . (split("/", $url))[-1]);

    if($resp =~ /4\d\d/){
        Carp::croak("You wanted to download $url but it doesn't exist or you don't have access! ($resp)");
    }
    if($resp =~/5\d\d/){
        Carp::croak("You wanted to download $url but there appears to be a problem with the host! ($resp)");
    }

    return RC_OK eq $resp;
}

sub open_file_for_writing {
    my ($self, $file) = @_;

    $self->validate_file_for_writing($file)
        or return;

    if (-e $file) {
        unless (unlink $file) {
            Carp::croak("Can't unlink $file: $!");
        }
    }

    return $self->_open_file($file, 'w');
}

sub open_gzip_file_for_reading {
    my ($self, $file) = @_;

    $self->validate_file_for_reading($file)
        or return;

    #check file type for gzip or symlink to a gzip
    my $file_type = $self->file_type($file);
    #debian bug #522441 - `file` can report gzip files as any of these....
    if ($file_type ne "gzip" && $file_type ne "Sun" && $file_type ne "Minix") {
        Carp::croak("File ($file) is not a gzip file");
    }

    my $pipe = "zcat ".$file." |";

    # _open_file throws its own exception if it doesn't work
    return $self->_open_file($pipe);
}

# Returns the file type, following any symlinks along the way to their target
sub file_type {
    my $self = shift;
    my $file = shift;

    $self->validate_file_for_reading($file);
    $file = $self->follow_symlink($file);

    my $result = `file -b $file`;
    my @answer = split /\s+/, $result;
    return $answer[0];
}

sub file_is_gzipped {
    my ($self, $filename) = @_;

    my $file_type = $self->file_type($filename);

    #NOTE: debian bug #522441 - `file` can report gzip files as any of these....
    if ($file_type eq "gzip" or $file_type eq "Sun" or $file_type eq "Minix") {
        return 1;
    } else {
        return 0;
    }
}

# Follows a symlink chain to reach the final file, accounting for relative symlinks along the way
sub follow_symlink {
    my $self = shift;
    my $file = shift;

    # Follow the chain of symlinks
    while (-l $file) {
        my $original_file = $file;
        $file = readlink($file);
        # If the symlink was relative, repair that
        unless (File::Spec->file_name_is_absolute($file)) {
            my $path = dirname($original_file);
            $file = join ("/", ($path, $file));
        }
        $self->validate_file_for_reading($file);
    }

    return $file;
}

sub get_file_extension_for_path {
    my $self = shift;
    my $path = shift;
    my ($extension) = $path =~ /(\.[^.]+)$/;
    return $extension;
}

#####
# Methods dealing with user names, groups, etc
#####
sub user_id {
    return $<;
}

sub username {
    my $class = shift;
    my $username = $ENV{'REMOTE_USER'} || getpwuid($class->user_id);
    return $username;
}

my $sudo_username = undef;
sub sudo_username {
    my $class = shift;
    unless(defined $sudo_username) {
        $sudo_username = $class->_sudo_username;
    }

    $sudo_username;
}

#split out for ease of testing
sub _sudo_username {
    my $class = shift;
    my $who_output = $class->cmd_output_who_dash_m || '';
    my $who_username = (split(/\s/,$who_output))[0] || '';

    my $sudo_username = $who_username eq $class->username ? '' : $who_username;
    $sudo_username ||= $ENV{'SUDO_USER'};

    return ($sudo_username || '');
}

sub current_user_is_admin {
    my $class = shift;
    return Genome::Sys->current_user_has_role('admin');
}

sub current_user_has_role {
    my ($class, $role_name) = @_;
    my $user = $class->current_user;
    return 0 unless $user;
    return $user->has_role_by_name($role_name);
}

sub current_user {
    my $class = shift;
    return Genome::Sys::User->get(username => $class->username);
}

sub cmd_output_who_dash_m {
    return `who -m`;
}

sub user_is_member_of_group {
    my ($class, $group_name) = @_;
    my $user = Genome::Sys->username;
    my $members = (getgrnam($group_name))[3];
    return ($members && $user && $members =~ /\b$user\b/);
}

#####
# Various utility methods
#####
sub open_browser {
    my ($class, @urls) = @_;
    for my $url (@urls) {
        if ($url !~ /:\/\//) {
            $url = 'http://' . $url;
        }
    }
    my $browser;
    if ($^O eq 'darwin') {
        $browser = "open";
    }
    elsif ($browser = `which firefox`) {

    }
    elsif ($browser = `which opera`) {

    }
    for my $url (@urls) {
        Genome::Sys->shellcmd(cmd => "$browser $url");
    }
    return 1;
}

sub shellcmd {
    # execute a shell command in a standard way instead of using system()\
    # verifies inputs and ouputs, and does detailed logging...

    # TODO: add IPC::Run's w/ timeout but w/o the io redirection...

    my ($self,%params) = @_;

    my %orig_params = %params;

    my $cmd                          = delete $params{cmd};
    my $output_files                 = delete $params{output_files};
    my $input_files                  = delete $params{input_files};
    my $output_directories           = delete $params{output_directories};
    my $input_directories            = delete $params{input_directories};
    my $allow_failed_exit_code       = delete $params{allow_failed_exit_code};
    my $allow_zero_size_output_files = delete $params{allow_zero_size_output_files};
    my $set_pipefail                 = delete $params{set_pipefail};
    my $allow_zero_size_input_files  = delete $params{allow_zero_size_input_files};
    my $skip_if_output_is_present    = delete $params{skip_if_output_is_present};
    my $redirect_stdout              = delete $params{redirect_stdout};
    my $redirect_stderr              = delete $params{redirect_stderr};
    my $dont_create_zero_size_files_for_missing_output =
        delete $params{dont_create_zero_size_files_for_missing_output};
    my $print_status_to_stderr       = delete $params{print_status_to_stderr};
    my $keep_dbh_connection_open     = delete $params{keep_dbh_connection_open};

    $set_pipefail = 1 if not defined $set_pipefail;
    $print_status_to_stderr = 1 if not defined $print_status_to_stderr;
    $skip_if_output_is_present = 1 if not defined $skip_if_output_is_present;
    if (%params) {
        my @crap = %params;
        Carp::confess("Unknown params passed to shellcmd: @crap");
    }

    my ($t1,$t2,$elapsed);

    # Go ahead and print the status message if the cmd is shortcutting
    if ($output_files and @$output_files) {
        my @found_outputs = grep { -e $_ } grep { not -p $_ } @$output_files;
        if ($skip_if_output_is_present
            and @$output_files == @found_outputs
        ) {
            $self->status_message(
                "SKIP RUN (output is present):     $cmd\n\t"
                . join("\n\t",@found_outputs)
            );
            return 1;
        }
    }

    my $old_status_cb = undef;
    unless  ($print_status_to_stderr) {
        $old_status_cb = Genome::Sys->message_callback('status');
        # This will avoid setting the callback to print to stderr
        # NOTE: we must set the callback to undef for the default behaviour(see below)
        Genome::Sys->message_callback('status',sub{});
    }

    if ($input_files and @$input_files) {
        my @missing_inputs;
        if ($allow_zero_size_input_files) {
            @missing_inputs = grep { not -e $_ } grep { not -p $_ } @$input_files;
        } else {
            @missing_inputs = grep { not -s $_ } grep { not -p $_ } @$input_files;
        }
        if (@missing_inputs) {
            Carp::croak("CANNOT RUN (missing input files):     $cmd\n\t"
                         . join("\n\t", map { -e $_ ? "(empty) $_" : $_ } @missing_inputs));
        }
    }

    if ($input_directories and @$input_directories) {
        my @missing_inputs = grep { not -d $_ } @$input_directories;
        if (@missing_inputs) {
            Carp::croak("CANNOT RUN (missing input directories):     $cmd\n\t"
                        . join("\n\t", @missing_inputs));
        }
    }

    # disconnect the db handle in case this is about to take awhile
    $self->disconnect_default_handles unless $keep_dbh_connection_open;

    if ($ENV{GENOME_SYS_PAUSE_SHELLCMD} and $cmd =~ $ENV{GENOME_SYS_PAUSE_SHELLCMD}) {
        my $file = '/tmp/GENOME_SYS_PAUSE.' . $$;
        $self->warning_message("RUN MANUALLY (and remove $file afterward): $cmd");
        Genome::Sys->write_file($file,$cmd . "\n");
        my $msg_time = time;
        while (-e $file) {
            sleep 3;
            if (time - $msg_time > (15*60)) {
                $msg_time = time;
                $self->warning_message("...waiting for $file to be removed after running command");
            }
        }
        $self->warning_message("resuming execution presuming the command was run manually");
    }
    else {
        $self->status_message("RUN: $cmd");

        $t1 = time();
        my $system_retval;
        eval {
                my ($restore_stdout);
                if ($redirect_stdout) {
                    no warnings 'once'; # OLDOUT is used only once, not
                    open(OLDOUT, '>&STDOUT') || die "Can't dup STDOUT: $!";
                    open(STDOUT, '>', $redirect_stdout) || die "Can't redirect stdout to $redirect_stdout: $!";
                    $restore_stdout = UR::Util::on_destroy(sub {
                        open(STDOUT, '>&OLDOUT');
                    });
                }

                my ($restore_stderr);
                if ($redirect_stderr) {
                    no warnings 'once'; # OLDERR is used only once, not
                    open(OLDERR, '>&STDERR') || die "Can't dup STDERR: $!";
                    open(STDERR, '>', $redirect_stderr) || die "Can't redirect stderr to $redirect_stderr: $!";
                    $restore_stderr = UR::Util::on_destroy(sub {
                        open(STDERR, '>&OLDERR');
                    });
                }

                # Set -o pipefail ensures the command will fail if it contains pipes and intermediate pipes fail.
                # Export SHELLOPTS ensures that if there are nested "bash -c"'s, each will inherit pipefail
                my $shellopts_part = 'export SHELLOPTS;';
                if ($set_pipefail) {
                    $shellopts_part = "set -o pipefail; $shellopts_part";
                } else {
                    $shellopts_part = "set +o pipefail; $shellopts_part";
                }

                {   # POE sets a handler to ignore SIG{PIPE}, that makes the
                    # pipefail option useless.
                    local $SIG{PIPE} = 'DEFAULT';
                    $system_retval = system('bash', '-c', "$shellopts_part $cmd");
                }

                print STDOUT "\n"; # add a new line so that bad programs don't break TAP, etc.
        };
        my $exception = $@;
        if ($exception) {
            Carp::croak("EXCEPTION RUNNING COMMAND. Failed to execute: $cmd\n\tException was: $exception");
        }
        my $child_exit_code = $system_retval >> 8;
        $t2 = time();
        $elapsed = $t2-$t1;

        if ( $system_retval == -1 ) {
            Carp::croak("ERROR RUNNING COMMAND. Failed to execute: $cmd\n\tError was: $!");

        } elsif ( $system_retval & 127 ) {
            my $signal = $system_retval & 127;
            my $withcore = ( $system_retval & 128 ) ? 'with' : 'without';
            Carp::croak("COMMAND KILLED. Signal $signal, $withcore coredump: $cmd");

        } elsif ($child_exit_code != 0) {
            if ($child_exit_code == 141) {
                my ($package, $filename, $line) = caller(0);
                my $msg = "SIGPIPE was recieved by command but IGNORED! cmd: '$cmd' in $package at $filename line $line";
                $self->error_message($msg);
            } elsif ($allow_failed_exit_code) {
                Carp::carp("TOLERATING Exit code $child_exit_code from: $cmd");
            } else {
                Carp::croak("ERROR RUNNING COMMAND.  Exit code $child_exit_code from: $cmd\nSee the command's captured STDERR (if it exists) for more information");
            }
        }
    }


    my @missing_output_files;
    if ($output_files and @$output_files) {
        @missing_output_files = grep { not -s $_ }  grep { not -p $_ } @$output_files;
    }
    if (@missing_output_files) {
        if ($allow_zero_size_output_files
            #and @$output_files == @missing_output_files
            # XXX This causes the command to fail if only a few of many files are empty, despite
            # that the option 'allow_zero_size_output_files' was given. New behavior is to warn
            # in either circumstance, and to warn that old behavior is no longer present in cases
            # where the command would've failed
        ) {
            if (@$output_files == @missing_output_files) {
                Carp::carp("ALL output files were empty for command: $cmd");
            } else {
                Carp::carp("SOME (but not all) output files were empty for command " .
                    "(PLEASE NOTE that earlier versions of Genome::Sys->shellcmd " .
                    "would fail in this circumstance): $cmd");
            }
            if ($dont_create_zero_size_files_for_missing_output) {
                @missing_output_files = (); # reset the list of missing output files
                @missing_output_files =
                    grep { not -e $_ }  grep { not -p $_ } @$output_files; # rescan for only missing files
            } else {
                for my $output_file (@missing_output_files) {
                    Carp::carp("ALLOWING zero size output file '$output_file' for command: $cmd");
                    my $fh = $self->open_file_for_writing($output_file);
                    unless ($fh) {
                        Carp::croak("failed to open $output_file for writing to replace missing output file: $!");
                    }
                    $fh->close;
                }
                @missing_output_files = ();
            }
        }
    }

    my @missing_output_directories;
    if ($output_directories and @$output_directories) {
        @missing_output_directories = grep { not -s $_ }  grep { not -p $_ } @$output_directories;
    }


    if (@missing_output_files or @missing_output_directories) {
        for (@$output_files) {
            if (-e $_) {
                unlink $_ or Carp::croak("Can't unlink $_: $!");
            }
        }
        Carp::croak("MISSING OUTPUTS! "
                    . join(', ', @missing_output_files)
                    . " "
                    . join(', ', @missing_output_directories));
    }
    unless  ($print_status_to_stderr) {
        # Setting to the original behaviour (or default)
        Genome::Sys->message_callback('status',$old_status_cb);
    }

    if ($ENV{GENOME_SYS_LOG_DETAIL}) {
        my $msg = encode_json({%orig_params, t1 => $t1, t2 => $t2, elapsed => $elapsed });
        Genome::Sys->debug_message(qq|$msg|)
    }

    return 1;

}

sub capture {
    my $class = shift;

    # lazy load so we don't break /gsc/bin/perl (until we have to)
    require IPC::System::Simple;
    return IPC::System::Simple::capture(@_);
}

sub disconnect_default_handles {
    my $class = shift;

    for my $ds (qw(Genome::DataSource::GMSchema)) {
        if($ds->has_default_handle) {
            $class->debug_message("Disconnecting $ds default handle.");
            $ds->disconnect_default_dbh();
        }
    }

    return 1;
}

sub retry {
    my %args = Params::Validate::validate(
        @_, {
            callback => { type => CODEREF },
            tries  => {
                type => SCALAR,
                callbacks => {
                    'is an integer' => sub { shift =~ /^\d+$/ },
                    'is greater than zero' => sub { shift > 0 },
                },
            },
            delay    => {
                type => SCALAR,
                callbacks => {
                    'is an integer' => sub { shift =~ /^\d+$/ },
                    'is greater than, or equal to, zero' => sub { shift >= 0 },
                },
            },
        },
    );

    my $rv;
    while ($args{tries} > 0) {
        $args{tries}--;
        $rv = $args{callback}->();
        last if $rv;
        sleep $args{delay};
    }

    return $rv;
}

1;

__END__

    methods => [
        dbpath => {
            takes => ['name','version'],
            uses => [],
            returns => 'FilesystemPath',
            doc => 'returns the path to a data set',
        },
        swpath => {
            takes => ['name','version'],
            uses => [],
            returns => 'FilesystemPath',
            doc => 'returns the path to an application installation',
        },
    ]

# until we get the above into ur...

=pod

=head1 NAME

Genome::Sys

=head1 VERSION

This document describes Genome::Sys version 0.9.0.1.

=head1 SYNOPSIS

 use Genome;
 my $dir = Genome::Sys->db_path('cosmic', 'latest');
 my $path = Genome::Sys->sw_path('htseq','0.5.3p9','htseq-count');

=head1 DESCRIPTION

Genome::Sys is a simple layer on top of OS-level concerns,
including those automatically handled by the analysis system,
like database cache locations.

=head1 METHODS

=head2 sw_path($pkg_name,$version) or sw_path($pkg_name,$version,$executable_basename)

Return the path to a given executable, library, or package.
The 3-parameter variation is only for packages which have multiple executables.

This is a wrapper for the OS-specific strategy for managing multiple versions of software packages,
(i.e. /etc/alternatives for Debian/Ubuntu)

The GENOME_SW environment variable contains a colon-separated lists of paths which this falls back to.
The default value is /var/lib/genome/sw/.

=head3 ex:

    $path = Genome::Sys->sw_path("tophat","2.0.7");

    $path = Genome::Sys->sw_path("htseq","0.5.3p9","htseq-count");

=head2 sw_versions($pkg_name) or sw_versions($pkg_name,$executable_basename);

Return a list of all installed versions of a given executable in order.

=head3 ex:

    @versions = Genome::Sys->sw_versions("tophat");

    @versions = Genome::Sys->sw_versions("htseq","htseq-count");

=head2 sw_version_path_map($pkg_name) or sw_version_path_map($pkg_name,$executable_basename)

Return a map of version numbers to executable paths.

=head3 ex:

    %map = Genome::Sys->sw_version_path_map("tophat");

    %map = Genome::Sys->sw_version_path_map("htseq","htseq-count");


=head2 db_path($name,$version)

Return the path to the preprocessed copy of the specified database.
(This is in lieu of a consistent API for the database in question.)

The GENOME_DB environment variable contains a colon-separated lists of paths which this falls back to.
The default value is /var/lib/genome/db/.

=head3 ex:

    my $dir1 = Genome::Sys->db_path('cosmic', 'latest');

=cut
