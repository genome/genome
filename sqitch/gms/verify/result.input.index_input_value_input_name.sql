-- Verify result.input.index_input_value_input_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'sri_iniv';

ROLLBACK;
