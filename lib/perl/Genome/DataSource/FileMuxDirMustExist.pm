package Genome::DataSource::FileMuxDirMustExist;

use Genome;
use Carp;
use IO::Dir;
use File::Basename;

class Genome::DataSource::FileMuxDirMustExist {
    is_abstract => 1,
};

my %verified_directories;
sub directory_must_exist {
    my $class = shift;
    my $path = shift;

    my $dir = File::Basename::dirname($path);
    unless (exists $verified_directories{$dir}) {
        $verified_directories{$dir} = IO::Dir->new($dir) ? 1 : 0;
    }
    unless ($verified_directories{$dir}) {
        my $type = UR::Object::Type->is_loaded(data_source_id => $class)
                    ||
                   __PACKAGE__->__meta__;
        Carp::croak('Containing directory does not exists for '
                    . $type->class_name . ": $path");
    }
    return $path;
}

1;
