-- Verify result.software_result.index_class_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'sr_cname';

ROLLBACK;
