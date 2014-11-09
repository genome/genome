-- Verify process.process.times_index

BEGIN;

SELECT 1/count(*) FROM pg_class
    WHERE relkind = 'i' and relname = 'process_process_created_at_idx';
SELECT 1/count(*) FROM pg_class
    WHERE relkind = 'i' and relname = 'process_process_started_at_idx';
SELECT 1/count(*) FROM pg_class
    WHERE relkind = 'i' and relname = 'process_process_ended_at_idx';

ROLLBACK;
