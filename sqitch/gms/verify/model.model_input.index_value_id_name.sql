-- Verify model.model_input.index_value_id_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'model_input_value_id_name_index';

ROLLBACK;
