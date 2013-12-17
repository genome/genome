-- Verify result_param

BEGIN;

SELECT software_result_id, param_name, param_value, name, value_class_name, value_id
FROM result.param
WHERE FALSE;

ROLLBACK;
