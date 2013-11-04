-- Revert instrument.fragment_library.index_sample_id

BEGIN;

DROP INDEX instrument.fragment_library_sample_id_index;

COMMIT;
