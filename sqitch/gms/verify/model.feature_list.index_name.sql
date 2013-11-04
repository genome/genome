-- Verify model.feature_list.index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'feature_list_name_index';

ROLLBACK;
