-- Verify disk.allocation.index_kilobytes_requested

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'kilobytes_requested_index';

ROLLBACK;
