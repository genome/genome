package Genome::Db;
use strict;
use warnings;
use Genome;

class Genome::Db {
    #is_abstract => 1,
    #subclassify_by => 'subclass_name',
    table_name => 'NA',
    id_by => ['id'],
    has_constant => [
        source_name         => { is => 'Text' },
        database_name       => { is => 'Text' },
        external_version    => { is => 'Text' },
        import_iteration    => { is => 'Text' },
        source_directory          => { is => 'FilesystemPath' },
        data_directory      => { is => 'FilesystemPath' },
        _subclass_name      => { is => 'ClassName' },
    ],
    doc => 'a versioned database' 
};

sub __extend_namespace__ {
    my ($this_class,$ext) = @_;
    return if $ext =~ /::/;
    my $new_class = $this_class . '::' . $ext;
    my $lc_ext = lc($ext);
    #print "Make $new_class\n";
    class {$new_class} {
        is => 'Genome::Db'
    };
}

sub __display_name__ {
    my $self = shift;
    return "DB(" . $self->id . ")";
}

sub __load__ {
    my ($class, $bx) = @_;

    my @rows;
    my @dirs = split(':',$ENV{GENOME_DB});
    for my $dir (@dirs) {
        my @db_dirs = glob("$dir/*");
        while (my $db_dir = shift @db_dirs) {
            unless (-d $db_dir) {
                #warn "not a directory db: $db_dir\n";
                next;
            }
            my $db = substr($db_dir,length($dir) + 1);
            
            my @version_dirs = glob("$db_dir/*");
            unless (grep { /\/latest$/ } @version_dirs) {
                #warn "no versions, retrying: @version_dirs\n";
                unshift @db_dirs, @version_dirs;
                next;
            }
            
            my ($source_name,$database_name);
            if ($db =~ /^(.*?)\/(.*)$/) {
                $source_name = $1;
                $database_name = $2;
            }
            else {
                $source_name = $db;
                $database_name = undef;
            }

            my $subclass_name = $source_name;
            $subclass_name =~ s/-/ /g;
            $subclass_name = Genome::Utility::Text::string_to_camel_case($subclass_name);
            $subclass_name = 'Genome::Db::' . $subclass_name;

            # warn " new db $db\n";
            for my $version_dir (@version_dirs) {
                my $version = substr($version_dir,length($db_dir) + 1);
                if ($version eq 'latest') {
                    next;
                }
                my $external_version;
                my $import_iteration;
                if ($version =~ /^(.*)[_-]v(.+)$/) {
                    $external_version = $1;
                    $import_iteration = $2;
                }
                elsif ($version =~ /^(.*)[\.](\d+)$/) {
                    $external_version = $1;
                    $import_iteration = $2;
                }
                else {
                    $external_version = $version;
                    $import_iteration = 1;
                }
                next unless -d $version_dir;
                next unless -e "$version_dir/.git";
                # warn "  new version $version\n";
                my $id = $source_name;
                $id .= '/' . $database_name if $database_name;
                $id .= '/';
                $id .= $external_version;
                $id .= '.' . $import_iteration if $import_iteration;
                push @rows, [$id,$source_name,$database_name,$external_version,$import_iteration,$dir,$version_dir,$subclass_name];
            }
        }
    }

    return ['id','source_name','database_name','external_version','import_iteration','source_directory','data_directory','_subclass_name'],\@rows;
}

1;

