-- Verify model_build_input

BEGIN;

SELECT build_id, value_class_name, value_id, name, filter_desc
FROM model.build_input
WHERE FALSE;

ROLLBACK;
