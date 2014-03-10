-- Verify disk.allocation.index_kilobytes_used

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'kilobytes_used_index';

ROLLBACK;
