-- Verify model_model_group

BEGIN;

SELECT id, name, user_name, uuid
FROM model.model_group
WHERE FALSE;

ROLLBACK;
