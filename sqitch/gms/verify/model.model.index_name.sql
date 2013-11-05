-- Verify model.model.index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'model_name_index';

ROLLBACK;
