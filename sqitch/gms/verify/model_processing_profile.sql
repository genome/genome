-- Verify model_processing_profile

BEGIN;

SELECT id, type_name, name, subclass_name
FROM model.processing_profile
WHERE FALSE;

ROLLBACK;
