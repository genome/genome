#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Model::MetagenomicComposition16s::Test;
use Test::More;

use_ok('Genome::Model::Event::Build::MetagenomicComposition16s::Reports');

my $model = Genome::Model::MetagenomicComposition16s::Test->model_for_sanger;
ok($model, 'got mc16s sanger model');
my $build = Genome::Model::Build->create( 
    id => -3388, 
    model => $model
);
ok($build, 'created build');
ok($build->get_or_create_data_directory, 'resolved data directory for build');
$build->create_subdirectories;

my $example_build = Genome::Model::MetagenomicComposition16s::Test->example_build_for_model($model);
ok($example_build, 'got example build');
ok(_link_example_data($build, $example_build), 'linked example data');

# set some values
$build->amplicons_attempted(5);
is($build->amplicons_attempted, 5, 'amplicons attempted set');
$build->amplicons_processed(5);
is($build->amplicons_processed, 5, 'amplicons processed set');
$build->amplicons_processed_success(1);
is($build->amplicons_processed_success, 1, 'amplicons processed success set');

$build->amplicons_classified(4);
is($build->amplicons_classified, 4, 'amplicons classified set');
$build->amplicons_classified_success(4);
is($build->amplicons_classified_success, 4, 'amplicons classified success set');

# run
my $reports = Genome::Model::Event::Build::MetagenomicComposition16s::Reports->create(build => $build, model => $model);
ok($reports, 'create');
$reports->dump_status_messages(1);
ok($reports->execute, 'execute');

# verify
my @reports = glob($build->reports_directory.'/*');
is(@reports, 2, "Created 2 reports");
ok(-s $build->reports_directory.'/Summary/report.xml', 'Created summary report');
ok(-s $build->reports_directory.'/Summary/report.html', 'Created summary report html');
ok(-s $build->reports_directory.'/Composition/report.xml', 'Created composition report');

#print $self->_build->data_directory."\n";<STDIN>;
done_testing();


#####################################

sub _link_example_data {
    my ($build, $example_build) = @_;

    for my $dir_to_link (qw/ chromat_dir edit_dir classification_dir /) {
        my $dest_dir = $build->$dir_to_link;
        Genome::Sys->validate_existing_directory($dest_dir) or die;
        my $source_dir = $example_build->$dir_to_link;
        my $dh = Genome::Sys->open_directory($source_dir) or die;

        $dh->read; $dh->read; # . and .. dirs
        while ( my $file = $dh->read ) {
            my $target = "$source_dir/$file";
            next if -d $target;
            my $link =  $dest_dir.'/'.$file;
            unless ( symlink($target, $link) ) {
                die "Can't symlink ($target) to ($link): $!.";
            }
        }
    }

    return 1;
}
