-- Verify result.software_result.index_lookup_hash_class_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'lookup_hash_class_name_idx';

ROLLBACK;
