-- Verify result.input.index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'result_input_name_index';

ROLLBACK;
