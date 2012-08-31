package Genome::Model::Tools::Assembly::MaizeBlast;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Assembly::MaizeBlast {
    is => 'Command',
    has => [
	string => {
	    is => 'Text',
	    doc => 'String of fasta to blast',
	},
	output_file => {
	    is => 'Text',
	    doc => 'User defined blast output file name',
	    is_optional => 1,
	},
    ],
};

sub help_brief {
    'Tools to blast string of fasta against maize db',
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
gmt assembly maize-blast --sequence -s agcgtgattccattggcattctagagtgaagtgatcgattc
gmt assembly maize-blast --sequence -s agcgtgattccattggcattctagagtgaagtgatcgattc --outut-file my_blast_out_file
EOS
}

sub execute {
    my $self = shift;

    my $db = $self->_validate_blast_db();
    unless ($db) {
	$self->error_message("Failed to validate maize bes db");
	return;
    }

    #create target file
    unlink 'blast_sequence.fasta';
    my $fh = Genome::Sys->open_file_for_writing('blast_sequence.fasta') ||
	return;
    $fh->print(">fasta\n".$self->string."\n");
    $fh->close;

    #run blast
    my $out = ($self->output_file) ? $self->output_file : Genome::Sys->username.$$;
    my $ec = system("blastn $db blast_sequence.fasta -V 100 -sort_by_highscore > $out");

    #just printing contents of blast output file
    my $fh2 = Genome::Sys->open_file_for_reading($out) ||
	return;
    my @contents = $fh2->getlines;
    $fh2->close;
    print map {$_} @contents;

    return 1;
}

sub _validate_blast_db {
    my $self = shift;

    my $db_dir = '/gscmnt/sata910/assembly/nthane/Maize_database';
    foreach ('xnd', 'xni', 'xns', 'xnt') {
	unless (-s $db_dir."/maize_454_db.$_") {
	    $self->error_message("Failed to find all components of blast db".
				 "\n\tMissing $db_dir/maize_454_db.$_");
	    return;
	}
    }
    return $db_dir.'/maize_454_db';
}

1;
