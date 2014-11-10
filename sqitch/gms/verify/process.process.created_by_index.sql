-- Verify process.process.created_by_index

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'process_created_by_idx';

ROLLBACK;
