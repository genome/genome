package Genome::Model::Tools::Db::ImportAlleleFreq;

use strict;
use warnings;

use Genome;                         # >above< ensures YOUR copy is used during development

class Genome::Model::Tools::Db::ImportAlleleFreq {
    is => 'Command',
    has => [
        group   => { is => 'String', doc => 'the name of the variant group to use' },
        dir     => { is => 'Directory', doc => 'the path to the directory of downloaded data.' },
    ],
};

sub sub_command_sort_position { 14 }

sub help_brief {
    "import variant allele frequencies"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt import-allele-freq --dir download-dir --group euro
EOS
}

sub help_detail {                           
    return <<EOS 
This imports HapMap data about variant allele frequencies.
EOS
}

sub execute {
    $DB::single = $DB::stopper;
    my $self = shift;
    my $dir = $self->dir;
    my $group_name = $self->group;

    unless (-d $dir) {
        $self->error_message("No directory $dir found!");
        return;
    }

    eval qq| use MPSampleData::VariationGroup; |;
    die $@ if $@;

    my @groups = MPSampleData::VariationGroup->search(
        group_name => $group_name        
    );
    if (@groups) {
        $self->error_message("Groups already present with name '$group_name': @groups");
        #return;
    }
    else {
        #(@groups) = MPSampleData::VariationGroup->insert({ group_name => $group_name });
    }

    my $expected_header = 'rs# chrom pos strand build center protLSID assayLSID panelLSID QC_code refallele refallele_freq refallele_count otherallele otherallele_freq otherallele_count totalcount';
    my @cols_in = split(/\s+/,$expected_header);
    $cols_in[0] =~ s/\#//g;

    my @cols_out = qw/chrom pos strand refallele refallele_freq refallele_count otherallele otherallele_freq otherallele_count/;

    my $sql = qq/
        select variation_id, start_, allele_string 
        from mg.variation v 
        join mg.chromosome c on c.chrom_id = v.chrom_id 
        where chromosome_name = ? 
        order by start_
    /;
    my $dbh = Genome::DataSource::GMSchema->get_default_dbh();
    $dbh->{RaiseError} = 1;
    my $sth = $dbh->prepare($sql);

    my @files = glob("$dir/allele_freqs*.gz");  
    for my $file (@files) {
        $self->status_message("checking $file");
        
        my ($expected_chrom) = ($file =~ /chr(.*?)_/); 
        unless ($expected_chrom) {
            $self->error_message("Failed to parse chromosome from file name $file!");
            return;
        }       

        my $fh = IO::File->new("zcat $file 2>&1 |");
        
        my $header = $fh->getline;
        chomp $header;
        unless ($header eq $expected_header) {
            die "File $file got\n'$header'\nexpected\n'$expected_header'\n";
        }
        
        $sth->execute($expected_chrom);

        my $cnt; 
        my $next_db_variant = $sth->fetchrow_hashref();
        ROW:
        while (my $row = $fh->getline) { 
            $cnt++;

            chomp $row;
            my %data;
            @data{@cols_in} = split(/\s/,$row);

            die "got $data{chrom} expected >$expected_chrom<" unless $data{chrom} eq 'chr' . $expected_chrom;

            my $allele_string = uc($data{refallele} . '/' . $data{otherallele});
            my $allele_freq = $data{otherallele_freq};
            my $pos = $data{pos};

            until ($next_db_variant->{START_} == $pos and $allele_string eq uc($next_db_variant->{ALLELE_STRING})) {
                if ($next_db_variant->{START_} > $pos) {
                    print "no variant at $pos for $expected_chrom and allele $allele_string\n";
                    next ROW;
                }
                $next_db_variant = $sth->fetchrow_hashref();
            }

            print join("\t", @data{@cols_out}, $expected_chrom, $pos, ':', $next_db_variant->{VARIATION_ID}, $allele_string, $allele_freq),"\n";
            last if $cnt==10;
        }

        $fh->close;
    }

    return 1;
}

1;

