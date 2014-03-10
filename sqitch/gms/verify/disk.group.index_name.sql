-- Verify disk.group.index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'group_name_index';

ROLLBACK;
