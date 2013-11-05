-- Verify result.input.index_input_value_input_name_software_result_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'sri_inivsri';

ROLLBACK;
