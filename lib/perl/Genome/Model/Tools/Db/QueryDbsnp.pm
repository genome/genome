package Genome::Model::Tools::Db::QueryDbsnp;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Db::QueryDbsnp {
    is => 'Command',                    
    has => [                                # specify the command's properties (parameters) <--- 
        chr     => { type => 'String',      doc => "give the chromosome name ie; 3 7 X" },
        coord     => { type => 'String',      doc => "give the NCBI Build 36 genomic coordinate" },
	query_result     => { type => 'String'      ,doc => "the query result not intended as an input", is_optional => 1  },
        dbsnp_version     => { type => 'String'      ,doc => "provide the dbsnp version ie 128 129 130" },
    ], 
};

sub help_brief {                            # keep this to just a few words <---
    "provide the chromosome and NCBI Build 36 coordinate and get dbsnp info in return"                 
}

sub help_synopsis {                         # replace the text below with real examples <---
    return <<EOS
gmt db query-dbsnp --chr=7 --coord=106311925 --dbsnp-version 130
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

	please provide the chromosome and NCBI Build 36 coordinate 
	for examplegt db query-dbsnp --chr=7 --coord=106311925 --dbsnp-version 130
	should result in rs2230460,'C/T',1 
	where the dbsnp is rs2230460 its alleles are C/T and 1 means its validated 
	0 inplace of 1 would indicate that it was not validated.

EOS
}

sub execute {                               # replace with real execution logic.
    my $self = shift;
    my $chr = $self->chr;
    my $coord = $self->coord;
    my $version = $self->dbsnp_version;
    unless ($chr && $coord && $version) { sub help_detail; }

    print "Running query-dbsnp $version command:\n";

    my $dbh = Genome::DataSource::Dwrac->get_default_handle;
    my $chrom_id = $dbh->prepare(qw/
				select seq_id from sequence_item si
				where sequence_item_type = 'chromosome sequence'
				and sequence_item_name = ?
				/);

    $chrom_id->execute('NCBI-human-build36-chrom' . $chr);
    
    
    my ($seq_id) = $chrom_id->fetchrow_array;
    
#---  statement handles allow you to deal with your results incrementally
    my $variation_exists = $dbh->prepare(qq/
					select ref_id,allele_description,is_validated from variation_sequence_tag vs
					join sequence_item si on si.seq_id=vs.vstag_id
					join sequence_tag st on st.stag_id = vs.vstag_id
					join sequence_correspondence scr on scr.scrr_id = vs.vstag_id
					join sequence_collaborator sc on sc.seq_id = vs.vstag_id
					where sc.collaborator_name = 'dbSNP'
					and sc.role_detail = '?'
					and seq2_start = ?
					and seq2_id = ? 
					/);
    
    $variation_exists->execute($version, $coord, $seq_id);
    
    my ($rs_id,$allele_description,$is_validated) = $variation_exists->fetchrow_array;
    my $result;
    if ($rs_id) {
	$result = qq($rs_id,$allele_description,$is_validated);
    } else {
	$result = 0;
	print qq(no dbSNP-$version found on NCBI-human-build36-chrom $chr at $coord\n);
    }
    $self->query_result($result);
    print qq($result\n);
    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;

