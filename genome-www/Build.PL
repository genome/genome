#!/usr/bin/env genome-perl
# Use local perl, not some perl on an application server!

use Config;
use Module::Build;

BEGIN {
  unshift @INC, "$ENV{PWD}/blib/lib";
}

my $class = Module::Build->subclass (
    class => 'Pod::Builder',
    code => <<'EOS',

      sub ACTION_clean {
        # FIXME: is this safe?
        use File::Path qw/rmtree/;
        rmtree "./_build";
        rmtree "./blib";
        rmtree "./gmt-bindoc";
        rmtree "./genome-bindoc";
        unlink "./Build";
        unlink "./MYMETA.yml";
      }

      sub ACTION_ur_docs {
        use File::Copy qw/copy/;
        $ENV{ANSI_COLORS_DISABLED} = 1;
        eval {
          local @INC = @INC;
          unshift @INC, 'blib/lib';
          die $@ if $@;
          eval "use Genome::Search";
          die $@ if $@;

          foreach my $exec ('genome','gmt') {
            UR::Namespace::Command::Update::Pod->execute(
                class_name => 'Genome::Search',
                executable_name => $exec,
                targets => ['Genome::Search'],
                output_path => 'cmd-bindoc'
            );
          }
        };
        die "failed to extract pod: $!: $@" if ($@);
      }

      sub ACTION_docs {
        my $self = shift;
        $self->depends_on('code');
        #$self->depends_on('ur_docs');
        $self->depends_on('manpages', 'html');
      }
       
      sub man1page_name {
        my ($self, $file) = @_;
        $file =~ s/.pod$//;
        return $self->SUPER::man1page_name($file);
      }

EOS
);

my $build = $class->new(
  module_name => 'Genome::Search',
  dist_version => 0.06,
  license => 'lgpl',
  requires => {
    'perl' => 'v5.8.7',
    'UR' => '0.29',
    'Carp' => '',
    'File::Basename' => '',
    'File::Temp' => '',
    'File::Find' => '',
    'IO::File' => '',
    'IO::String' => '',
    'IO::Socket' => '',
    'Sys::Hostname' => '',
    'Plack' => '',
    'Cache::Memcached' => '',
    'AnyEvent' => '',
    'Net::Server' => '',
    'Starman' => '',
    'Path::Class' => '',
    'XML::LibXML' => '',
    'XML::LibXSLT' => '',
    'XML::Simple' => '',
  },
  cpan_client => 'cpanm',
  etc_files => {
      'etc/solr/solr-tomcat.xml' => 'solr/solr-tomcat.xml',
      'etc/solr/tomcat.policy' => 'solr/tomcat.policy',
      'etc/solr/web.xml' => 'solr/web.xml',
      'etc/solr/conf/admin-extra.html' => 'solr/conf/admin-extra.html',
      'etc/solr/conf/elevate.xml' => 'solr/conf/elevate.xml',
      'etc/solr/conf/mapping-ISOLatin1Accent.txt' => 'solr/conf/mapping-ISOLatin1Accent.txt',
      'etc/solr/conf/protwords.txt' => 'solr/conf/protwords.txt',
      'etc/solr/conf/schema.xml' => 'solr/conf/schema.xml',
      'etc/solr/conf/scripts.conf' => 'solr/conf/scripts.conf',
      'etc/solr/conf/solrconfig.xml' => 'solr/conf/solrconfig.xml',
      'etc/solr/conf/spellings.txt' => 'solr/conf/spellings.txt',
      'etc/solr/conf/stopwords.txt' => 'solr/conf/stopwords.txt',
      'etc/solr/conf/synonyms.txt' => 'solr/conf/synonyms.txt',
      'etc/solr/conf/xslt/luke.xsl' => 'solr/conf/xslt/luke.xsl',
      'etc/genome_solr.xml' => 'tomcat/genome_solr.xml',
      'etc/genome' => 'apache/genome',
      'etc/init/genome-task-runner.conf' => 'init/genome-task-runner.conf',
      'etc/init/genome-search-index-queue.conf' => 'init/genome-search-index-queue.conf',
      'etc/cron.daily/genome-task-cleanup.sh' => 'cron.daily/genome-task-cleanup.sh',
  },
  install_path => {
      'solr' => '/var/cache/genome/solr/',
      'apache' => '/etc/apache2/sites-available/',
      'tomcat' => '/etc/tomcat6/Catalina/localhost/',
      'init' => '/etc/init',
      'cron.daily' => '/etc/cron.daily',
  },
);

$build->add_build_element('psql');
$build->add_build_element('etc');
$build->add_build_element('solr');
$build->add_build_element('tomcat');
$build->add_build_element('apache');
$build->add_build_element('psgi');
$build->add_build_element('css');
$build->add_build_element('txt');
$build->add_build_element('png');
$build->add_build_element('html');
$build->add_build_element('jpg');
$build->add_build_element('gif');
$build->add_build_element('js');
$build->add_build_element('as');
$build->add_build_element('swf');
$build->add_build_element('sql');
$build->add_build_element('php');
$build->add_build_element('xsl');
$build->add_build_element('xml');
$build->add_build_element('conf');
$build->add_build_element('sh');
$build->create_build_script;
