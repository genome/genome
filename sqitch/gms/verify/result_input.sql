-- Verify result_input

BEGIN;

SELECT software_result_id, input_name, input_value, name, value_class_name, value_id
FROM result.input
WHERE FALSE;

ROLLBACK;
