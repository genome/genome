
# review jlolofie
# 1. last matching line should be the new min everytime (for each chromosome file)
# 2. this should be benchmarked against another intersect tool to figure out which is faster-
#    depends on the density of the data

package Genome::Model::Tools::Annotate::LookupVariants;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use IO::File;
use Carp;

use Bio::DB::Fasta;

class Genome::Model::Tools::Annotate::LookupVariants {
    is  => 'Genome::Model::Tools::Annotate',
    has => [
        variant_file => {
            type     => 'Text',
            is_input => 1,
            doc      =>
                "File of variants. TSV (sorted by chromosome,start): chromosome_name start stop reference variant",
        },
        output_file => {
            type      => 'Text',
            is_input  => 1,
            is_output => 1,
            default   => "STDOUT",
            doc       => "default is STDOUT",
        },
    ],
    has_optional => [
        _output_filehandle => {
            type      => 'SCALAR',
        },
        _last_data_line_number => {
            type      => 'SCALAR',
            doc => 'The number of lines in the input',
        },
        filter_out_submitters => {
            type     => 'Text',
            is_input => 1,
            doc      =>
                'Comma separated (no spaces allowed) list of submitters to IGNORE from dbsnp',
        },
        dbSNP_path => {
            type     => 'Text',
            is_optional => 1,
            doc      => "path to dbSNP files broken into chromosome",
            example_values => ['/gscmnt/sata835/info/medseq/model_data/2857166586/ImportedVariations/tmp'],
        },
        dbSNP_version => {
            type    => 'Int',
            is_input => 1,
            is_optional => 1,
            default => 130,
            doc     => 'Version of dbSNP to use. Defaults to 130.',
        },
        report_mode => {
            type     => 'Text',
            is_input => 1,
            default  => 'novel-only',
            doc      =>
           '     - novel-only (DEFAULT VALUE) prints lines from variant_file that are not found in dbSNP
     - known-only => prints lines from variant file that are found in dbSNP
     - full => gives back each line from your variant file with matched lines from dbsnp 
     - gff => gives back all the lines that match the range of coordinates in your variant file',
        },
        index_fixed_width => {
            type     => 'Int',
            calculate_from  => ['dbSNP_version'],
            calculate      => sub { # if we want to use a more efficient index file, this will be easy to change.
                                    return 10;
                                },
        },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
        },
        append_rs_id => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'append rs_id from dbSNP at end of each matching row'
        },
        append_allele=> {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'append allele from dbSNP to end of output'
        },
        append_population_allele_frequencies=> {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'use in junction with report-mode full or gff, to append population allele frequencies to end of output'
        },
        frequency_path => {
            type     => 'Text',
            is_optional => 1,
            doc      => "path to allele frequency directory... this defaults to human from dbsnp 129 data. Mouse is available at /gscmnt/sata835/info/medseq/model_data/2857225771/ImportedVariations/frequencies",
            example_values => ['/gscmnt/sata835/info/medseq/model_data/2857282699/ImportedVariations/frequencies'],
        },
        organism => {
            type  =>  'String',
            doc   =>  "provide the organism if annotating mouse",
            valid_values => ['human', 'mouse'],
            is_optional  => 1,
            default => 'human',
        },
        require_allele_match => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'filter out SNPs where allele doesnt match'
        }
    ],
};


sub help_synopsis { 
    return <<EOS
gmt annotate lookup-variants --variant-file snvs.csv --output-file novel_variants.csv
EOS
}

sub help_detail {
    return <<EOS
By default, takes in a file of variants and filters out variants that are already known to exist.
EOS
}

sub execute { 

    my ($self) = @_;

    my $dbsnp_dir= $self->dbSNP_path;
    unless($dbsnp_dir) {
        $self->error_message('Not given dbsnp path');
        Carp::confess $self->error_message;
    }
    $self->dbSNP_path($dbsnp_dir);

    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    my $variant_file = $self->variant_file;
    open(my $in, $variant_file) || die "cant open $variant_file";

    my $fh = $self->get_output_fh() || die 'no output filehandle';
    $self->_output_filehandle($fh);

    while (my $line = <$in>) {
        if(substr($line, 0, 5) ne "chrom") {
            $self->print_matches($line);
        }
    }

    close($in);
    $fh->close;
    
    for my $chr (keys %{ $self->{'filehandles'} }) {
        $self->{'filehandles'}->{$chr}->close;
    }

    return 1;
}

sub print_matches {

    my ($self, $line) = @_;
    my ($chr, $start, $stop, $allele1, $allele2) = split(/\t/,$line);

    my @matches = $self->find_all_matches($line);
    @matches = map { $self->filter_by_type($_) } @matches; # SNPs only

    if ($self->filter_out_submitters()) {
        @matches = map { $self->filter_by_submitters($_) } @matches;
    }

    if ($self->require_allele_match()) {
        chomp($allele2);
        @matches = map { $self->filter_by_allele($_, $allele1, $allele2) } @matches;
    }

    my $fh = $self->_output_filehandle() || die 'no output_filehandle';
    my $report_mode = $self->report_mode();
    if (($report_mode eq 'known-only')&&(@matches)) {

        if ($self->append_rs_id() || $self->append_allele() ) {

            my ($dbsnp_fh, $index) = $self->get_fh_for_chr($chr);
            my $snp = parse_dbsnp_line($matches[0]);

            if ($self->append_rs_id()) {
                my $rs_id = $snp->{'rs_id'};
                chomp($line);
                $line = sprintf("%s\t%s\n",$line,$rs_id); 
            }

            if ($self->append_allele()) {
                my $allele = $snp->{'ds_allele'};
                chomp($line);
                $line = sprintf("%s\t%s\n",$line,$allele); 
            }
        }

        $fh->print($line);
    } elsif (($report_mode eq 'novel-only')&& (scalar @matches == 0)) {
        $fh->print($line);

    } elsif ($report_mode eq 'gff') {

        my @matches;
        for my $n ($start..$stop) {
            my $new_line = qq($chr\t$n\t$n);
            @matches = $self->find_all_matches($new_line);
            for my $snp_line (@matches) {

                my (@rs_ids,@submitters,@matchs);
                my $snp = parse_dbsnp_line($snp_line);

                my $rs_id = $snp->{'rs_id'};
                my $allele = $snp->{'ds_allele'};
                my $submitter = $snp->{'ds_submitter'};

                my $ds_type = $snp->{'ds_type'};
                my $ds_start = $snp->{'ds_start'};
                my $ds_stop = $snp->{'ds_stop'};
                my $strain = $snp->{'strain'};
                my $is_validated = $snp->{'is_validated'};
                my $is_validated_by_allele = $snp->{'is_validated_by_allele'};
                my $is_validated_by_cluster = $snp->{'is_validated_by_cluster'};
                my $is_validated_by_frequency = $snp->{'is_validated_by_frequency'};
                my $is_validated_by_hap_map = $snp->{'is_validated_by_hap_map'};
                my $is_validated_by_other_pop = $snp->{'is_validated_by_other_pop'};

                my @validated;
                if ($is_validated_by_allele == 1) {
                    push(@validated,"by2Hit2Allele");
                }
                if ($is_validated_by_cluster == 1) {
                    push(@validated,"byCluster");
                }
                if ($is_validated_by_frequency == 1) {
                    push(@validated,"byFrequency");
                }
                if ($is_validated_by_hap_map == 1) {
                    push(@validated,"byHapMap");
                }
                if ($is_validated_by_other_pop == 1) {
                    push(@validated,"byOtherPop");
                }
                unless (@validated) {@validated = "not_validated";}
                my $validation = join ':' , @validated;



                my $gff_line = qq(Chromosome$chr\tdbsnp_130\t$ds_type\t$ds_start\t$ds_stop\t.\t$strain\t.\t$rs_id \; Alleles \"$allele\" ; validation \"$validation\" ; submitter \"$submitter\");



                if ($self->append_population_allele_frequencies) {
                    my $freq = &get_frequencies($self,$rs_id);
                    if ($freq) {
                        chomp $freq;
                        $gff_line = "$gff_line ; $freq";
                    }
                }

                $fh->print("$gff_line\n");
            }
        }

    } elsif ($report_mode eq 'full') {
        my $report_line;
        chomp($line);

        if (@matches) {
            my (@rs_ids,@submitters,@matchs,@validation);
            for my $snp_line (@matches) {
                my $snp = parse_dbsnp_line($snp_line);
                my $rs_id = $snp->{'rs_id'};
                my $allele = $snp->{'ds_allele'};
                my $submitter = $snp->{'ds_submitter'};
                my $is_validated = $snp->{'is_validated'};
                my $is_validated_by_allele = $snp->{'is_validated_by_allele'};
                my $is_validated_by_cluster = $snp->{'is_validated_by_cluster'};
                my $is_validated_by_frequency = $snp->{'is_validated_by_frequency'};
                my $is_validated_by_hap_map = $snp->{'is_validated_by_hap_map'};
                my $is_validated_by_other_pop = $snp->{'is_validated_by_other_pop'};

                #my $validation = "$is_validated\-$is_validated_by_allele\-$is_validated_by_cluster\-$is_validated_by_frequency\-$is_validated_by_hap_map\-$is_validated_by_other_pop";
                #chomp $validation;

                my @validated;
                if ($is_validated_by_allele == 1) {
                    push(@validated,"by2Hit2Allele");
                }
                if ($is_validated_by_cluster == 1) {
                    push(@validated,"byCluster");
                }
                if ($is_validated_by_frequency == 1) {
                    push(@validated,"byFrequency");
                }
                if ($is_validated_by_hap_map == 1) {
                    push(@validated,"byHapMap");
                }
                if ($is_validated_by_other_pop == 1) {
                    push(@validated,"byOtherPop");
                }
                unless (@validated) {@validated = "not_validated";}

                #my $validation = join ';' , @validated;
                #push(@validation,$validation) unless grep (/$validation/ , @validation);

                push(@rs_ids,$rs_id) unless grep (/$rs_id/ , @rs_ids);
                push(@submitters,$submitter) unless grep (/$submitter/ , @submitters);

                my $match;
                if ($allele1 && $allele2) {
                    chomp($allele2);
                    $match = $self->filter_by_allele($snp_line, $allele1, $allele2);
                    if ($match) {$match = "$allele:allele_match";} else {$match = "$allele:no_match";}
                } else {
                    $match = $allele;
                }
                push(@matchs,$match) unless grep (/$match/ , @matchs);

                if ($match =~ /no_match/) {push(@validated,"alternate_allele");}
                my $validation = join ';' , @validated;
                push(@validation,$validation) unless grep (/$validation/ , @validation);

            }
            my $rs_id = join ':' , @rs_ids;
            my $submitter = join ':' , @submitters;
            my $match = join ':' , @matchs;
            my $validation  = join ':' , @validation;

            my $frequencies;
            if ($self->append_population_allele_frequencies) {
                my @freq;
                for my $rsid (@rs_ids) {
                    my $freq = &get_frequencies($self,$rs_id);
                    if ($freq) {
                        push(@freq,$freq);
                    }
                }
                if (@freq) {
                    $frequencies = join ':' , @freq;
                } else {
                    $frequencies = "-";
                }
                $report_line = sprintf("%s\t%s\t%s\n",$line,"$rs_id,$submitter,$match,$validation",$frequencies);
            } else {
                $report_line = sprintf("%s\t%s\n",$line,"$rs_id,$submitter,$match,$validation");
            }

        } else {
            my $match = "no_hit";
            $report_line = sprintf("%s\t%s\n",$line,$match); 
        }

        $fh->print($report_line);
    }
}

sub get_output_fh {
    my $self = shift;

    my $output_file = $self->output_file();
    die 'no output_file!' if !$output_file;

    if ($output_file eq 'STDOUT') {
        return 'STDOUT';
    }

    my $fh = IO::File->new(">" . $output_file);
    return $fh;
}

sub filter_by_allele {

    my ($self, $line, $ref, $variant) = @_;

    my $snp = parse_dbsnp_line($line);    
    my $ds_allele = $snp->{'ds_allele'};
    #my ($dbsnp_allele1,$dbsnp_allele2) = split(/\//,$ds_allele);
    my @dbsnp_allele_array = split(/\//,$ds_allele);
    my $array_n = @dbsnp_allele_array;

    use Genome::Info::IUB;
   
    # if variant is not expanded, include ref in alpha order 
    my ($a1, $a2) = Genome::Info::IUB->variant_alleles_for_iub($ref,$variant); ##if $variant eq "N" there would be an a3

#TODO: finish work here

    my @vars;
    if ($a1) {push(@vars,$a1);}
    if ($a2) {push(@vars,$a2);}

    my ($rm,$vm);
    for my $var (@vars) {
        for my $n (1..$array_n) {
            my $m = $n - 1;
            my $dbsnp_allele = $dbsnp_allele_array[$m];
            if ($ref eq $dbsnp_allele) {
                $rm = $dbsnp_allele;
            } elsif ($var eq $dbsnp_allele) {
                $vm = $dbsnp_allele;
            }
        }
        unless ($rm && $vm) {
            undef($rm);
            undef($vm);
            for my $n (1..$array_n) {
                my $m = $n - 1;
                my $dbsnp_allele = $dbsnp_allele_array[$m];
                my $rev_dbsnp_allele = &reverse_complement ($dbsnp_allele); 
                if ($ref eq $rev_dbsnp_allele) {
                    $rm = $rev_dbsnp_allele;
                } elsif ($var eq $rev_dbsnp_allele) {
                    $vm = $rev_dbsnp_allele;
                }
            }
        }
    }

    if ($rm && $vm) {
        return $line;
    } else {
        return;
    }
}

sub reverse_complement {
    my ($sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/aAcCtTgG/tTgGaAcC/;
    return $sequence;
}

sub get_frequencies {

    my ($self,$rs_id) = @_;

    # Temporary fix for an immediate need... these files should be moved to a build directory for the models that used to hold them in model directories
    # These models were ImportedVariation models with the names dbSNP-human-129 and dbSNP-mouse-10090 and resided in $model->data_directory."/ImportedVariations/frequencies";
    my $RsDir = $self->frequency_path;

    my $rsdb = Bio::DB::Fasta->new($RsDir);
    chmod(0666, $RsDir . "/directory.index"); #change the index permissions so the next guy can rebuild it
    my $freq = $rsdb->seq($rs_id, 1 => 6000);

    return unless $freq;
    return $freq;

}

sub filter_by_submitters {

    my ($self, $line) = @_;

    # NOTE: submitters are 1 per line in dbsnp data source files,
    # but are comma separated list on command line

    my $snp = parse_dbsnp_line($line);    
    my $ds_submitter = $snp->{'ds_submitter'};

    my $submitters_str = $self->filter_out_submitters();
    return $line if !$submitters_str;


    my @filter_submitters = split(/,/,$submitters_str);
    if (grep /^$ds_submitter$/, @filter_submitters) {
        return;
    }

    return $line;
}


sub filter_by_type {

    my ($self, $line) = @_;
    # returns $line if its a SNP

    my $snp = parse_dbsnp_line($line);

    if ($snp->{'ds_type'} eq 'SNP'
        && $snp->{'ds_start'} == $snp->{'ds_stop'}) {
        return $line;
    }

    return;
}

sub find_all_matches {

    # TODO: the problem is we only return position, not chromosome, etc

    my ($self, $line) = @_;
    my @matches;

    my $pos = $self->find_a_matching_pos($line);
    if (defined ($pos)) {
        my ($chr, $start, $stop) = split(/\t/,$line);
        @matches = $self->find_matches_around($chr, $pos);
    }

    return @matches;
}

sub find_matches_around {

    my ($self, $chr, $pos) = @_;
    my $variant = {};
    my (@forward, @backward);

    return unless grep {$_ eq $chr } ( 1..22, 'X', 'Y');

    my ($fh, $index) = $self->get_fh_for_chr($chr);

    my $cur = $pos;
    my $original_line = $self->get_line($fh, $index, $cur);

    my $snp = parse_dbsnp_line($original_line);
    my $ds_start = $snp->{'ds_start'};
    my $start = $ds_start;

    push @forward, $original_line;
   
    # go forward 
    while ($cur == $pos || $start == $ds_start) {
        $cur++;
        last if ($cur > $self->_last_data_line_number);
        
        my $forward_line = $self->get_line($fh, $index, $cur);
        last if !$forward_line;

        my $forward_snp = parse_dbsnp_line($forward_line);
        $ds_start = $forward_snp->{'ds_start'};
        if ($start == $ds_start) {
            push @forward, $forward_line;
        }
    }

    # reset and go backwards
    $ds_start = $start;
    $cur = $pos; 
    while ($cur == $pos || $start == $ds_start) {
        $cur--;
        last if ($cur < 0);

        my $reverse_line = $self->get_line($fh, $index, $cur);
        last if !$reverse_line;

        my $reverse_snp = parse_dbsnp_line($reverse_line);
        $ds_start = $reverse_snp->{'ds_start'};

        if ($start == $ds_start) {
            push @backward, $reverse_line;
        }
    } 

    return (reverse @backward, @forward);
}

sub find_a_matching_pos {

    my ($self, $line) = @_;

    my ($chr, $start, $stop) = split(/\t/,$line);

    return unless grep {$_ eq $chr } ( 1..22, 'X', 'Y');

    my ($fh, $index) = $self->get_fh_for_chr($chr);
    my $match_count = 0;
    my $size = <$index>; chomp($size);
    my $min = 0;
    my $max = $size - 1;
    $self->_last_data_line_number($max);

    while($min <= $max) {

        my $cur += $min + int(($max - $min) / 2);

        my $line = $self->get_line($fh, $index, $cur);
        my $snp = parse_dbsnp_line($line);
        my $ds_start = $snp->{'ds_start'};

        if ($start > $ds_start) {
            $min = $cur + 1;
        } elsif ($start < $ds_start) { 
            $max = $cur - 1;
        } else {
            return $cur;
        }
    }

    return;
}

sub get_line {

    my ($self, $fh, $index, $line_number) = @_;
    my $fixed_width = $self->index_fixed_width();

    # add fixed width to account for index header
    my $index_pos = $line_number * $fixed_width + $fixed_width;
    seek($index, $index_pos, 0);
    my $pos = <$index>; chomp($pos);
   
    seek($fh, $pos, 0);
    my $line = <$fh>;
    return $line; 
}

sub get_fh_for_chr {

    my ($self, $chr) = @_;

    my $dbSNP_path = $self->dbSNP_path();

    my $organism = $self->organism;
    my $model;

    #TODO should this simply check to see if organism eq mouse?

    #if ($dbSNP_path eq "/gsc/var/lib/import/dbsnp/130/tmp/" && $organism eq "mouse") {
    if ($organism eq "mouse") {
        unless($model = Genome::Model::ImportedVariations->get(name => 'dbSNP-mouse-10090')){
            $self->error_message("Could not locate ImportedVariations model with dbSNP-mouse-10090 name");
            die $self->error_message;
        }
#$dbSNP_path = "/gsc/var/lib/import/dbsnp/mouse/10090/";
#        $dbSNP_path = $model->data_directory."/ImportedVariations/";
        die 'Error: please email $ENV{GENOME_EMAIL_PIPELINE} -- mouse dbsnp used to be stored in a model data directory';
        # This will probably work if you just change the dbsnp path? Would have to test it, the files look different
    }

    my ($fh, $index);

    if (!$self->{'filehandles'}->{$chr}) {
        my $dbSNP_filename = join('',  'variations_', $chr, '.csv');
        my $dbSNP_pathname = join('/',$dbSNP_path,$dbSNP_filename); 
        die "cant open dbSNP_pathname = $dbSNP_pathname" if ! -e $dbSNP_pathname ;

        my $index_filename = join('',  'variations_', $chr, '.csv.index');
        my $index_pathname = join('/',$dbSNP_path,$index_filename); 
        die "cant open index_pathname = $index_pathname" if ! -e $index_pathname ;

        open($fh, $dbSNP_pathname);
        $self->{'filehandles'}->{$chr} = $fh;

        open($index, $index_pathname);
        $self->{'index_filehandles'}->{$chr} = $index;
    } else {
        $fh = $self->{'filehandles'}->{$chr};
        seek($fh, 0, 0);

        $index = $self->{'index_filehandles'}->{$chr};
        seek($index, 0, 0);
    }

    die "no filehandle for $chr" if !$fh || !$index;

    return ($fh, $index);
}


sub parse_dbsnp_line {

    my ($line) = @_;

    my $snp = {};
    my @parts = split(/\t/,$line);

    my @keys = qw(
        ds_id
        ds_allele 
        ds_type
        ds_chr
        ds_start
        ds_stop
        ds_submitter
        rs_id
        strain
        is_validated
        is_validated_by_allele
        is_validated_by_cluster
        is_validated_by_frequency
        is_validated_by_hap_map
        is_validated_by_other_pop
    );

    my $i = 0;
    for my $key (@keys) {
        $snp->{$key} = $parts[$i];
        $i++;
    }

    return $snp;
}

1;






=pod

=head1 Name

Genome::Model::Tools::Annotate::LookupVariations

=head1 Synopsis

By default, takes in a file of variants and filters out variants that are already known to exist.

=head1 Usage

    $ gmt annotate lookup-variants --variant-file snvs.csv --output-file novel_variants.csv
 
=cut


