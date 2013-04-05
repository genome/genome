package Genome::Db::Base;
use strict;
use warnings;
use Genome;

class Genome::Db::Base {
    id_by => [
        database_name   => { is => 'Text' },
        version         => { is => 'Text' },
    ],
    has_constant => [
        source_dir      => { is => 'FilesystemPath' },
        dataset_dir     => { is => 'FilesystemPath' },
    ],
    composite_id_separator => ':',
    doc => 'a versioned database' 
};

sub __display_name__ {
    my $self = shift;
    return $self->class . " (" . $self->id . ")";
}

sub __load__ {
    my ($class, $bx) = @_;

    $DB::single = 1;

    my @dirs = split(':',$ENV{GENOME_DB});

    my $len_db = length($ENV{GENOME_DB})+1;
    my @db1 = map { glob("$_/*") } @dirs;

    my @rows;
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
            
            # warn " new db $db\n";
            for my $version_dir (@version_dirs) {
                my $version = substr($version_dir,length($db_dir) + 1);
                if ($version eq 'latest') {
                    next;
                }
                next unless -d $version_dir;
                next unless -e "$version_dir/.git";
                # warn "  new version $version\n";
                my $id = join('/',$db,$version);
                push @rows, [$id,$db,$version,$dir,$version_dir];
            }
        }
    }

    return ['id','database_name','version','source_dir','dataset_dir',],\@rows;
}

1;

