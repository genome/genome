package BAP::DB::SequenceSet;

use base 'BAP::DB::DBI';


__PACKAGE__->table('sequence_set');
__PACKAGE__->columns(
                     'All' => qw(
                                 sequence_set_id
                                 sequence_set_name
                                 organism_id
                                 software_version
                                 data_version
                                )
                    );
__PACKAGE__->sequence('sequence_set_id_seq');
__PACKAGE__->has_many('sequences' => 'BAP::DB::Sequence');

1;
