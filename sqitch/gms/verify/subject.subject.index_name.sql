-- Verify subject.subject.index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'subject_name_index';

ROLLBACK;
