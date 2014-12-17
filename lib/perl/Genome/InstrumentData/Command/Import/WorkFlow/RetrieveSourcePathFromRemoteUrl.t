#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromRemoteUrl') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v5') or die;
my $source_path = $test_dir.'/input.bam';
my $source_md5_path = $test_dir.'/input.bam.md5';

# classes to overload LWPUserAgent functionality
class Genome::LWPUserAgent {};
sub Genome::LWPUserAgent::get { 
    my %h; return bless \%h, 'Genome::LWPUserAgentResponse'; 
}

class Genome::LWPUserAgentResponse {};
sub Genome::LWPUserAgentResponse::is_success { return 1; };
sub Genome::LWPUserAgentResponse::headers { my %h; return bless \%h, 'Genome::LWPUserAgentHeaders'; };

class Genome::LWPUserAgentHeaders {};
my $cnt = 0;
sub Genome::LWPUserAgentHeaders::content_length { return ++$cnt == 1 ? 20037 : 121; };

Sub::Install::reinstall_sub({
        code => sub {
            my %h;
            return bless \%h, 'Genome::LWPUserAgent';
        },
        into => 'LWP::UserAgent',
        as => 'new',
    });


my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromRemoteUrl->execute(
    working_directory => $test_dir,
    source_path => $source_path,
);
ok($cmd->result, 'execute');
my $destination_path = $cmd->destination_path;
is($destination_path, $test_dir.'/input.bam', 'destination path named correctly');
ok(-s $destination_path, 'destination path exists');
my $original_md5_path = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->original_md5_path_for($destination_path);
ok(-s $original_md5_path, 'destination md5 exists');

done_testing();
