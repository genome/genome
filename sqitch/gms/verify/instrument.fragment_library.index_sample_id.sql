-- Verify instrument.fragment_library.index_sample_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'fragment_library_sample_id_index';

ROLLBACK;
