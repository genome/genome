-- Verify process.process.status_index

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'process_status_idx';

ROLLBACK;
