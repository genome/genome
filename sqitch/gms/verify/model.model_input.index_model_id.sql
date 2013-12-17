-- Verify model.model_input.index_model_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'model_input_model_id_index';

ROLLBACK;
