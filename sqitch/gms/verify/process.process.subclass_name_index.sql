-- Verify process.process.subclass_name_index

BEGIN;

SELECT 1/count(*) FROM pg_class
    WHERE relkind = 'i' and relname = 'process_subclass_name_idx';

ROLLBACK;
