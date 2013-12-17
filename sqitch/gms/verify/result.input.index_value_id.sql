-- Verify result.input.index_value_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'result_input_value_id_index';

ROLLBACK;
