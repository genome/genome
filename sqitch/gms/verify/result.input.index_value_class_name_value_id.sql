-- Verify result.input.index_value_class_name_value_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'result_input_value_class_id_index';

ROLLBACK;
