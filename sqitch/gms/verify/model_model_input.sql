-- Verify model_model_input

BEGIN;

SELECT model_id, value_class_name,value_id, name,filter_desc
FROM model.model_input
WHERE FALSE;

ROLLBACK;
