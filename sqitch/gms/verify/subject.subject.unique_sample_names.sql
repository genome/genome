-- Verify subject.subject.unique_sample_names

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' AND relname = 'unique_sample_name_index';

ROLLBACK;
