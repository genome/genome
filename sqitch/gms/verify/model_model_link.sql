-- Verify model_model_link

BEGIN;

SELECT to_model_id, from_model_id, role
FROM model.model_link
WHERE FALSE;

ROLLBACK;
