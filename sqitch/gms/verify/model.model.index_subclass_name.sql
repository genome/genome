-- Verify model.model.index_subclass_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'model_subclass_index';

ROLLBACK;
