-- Verify model_build_link

BEGIN;

SELECT to_build_id, from_build_id, role
FROM model.build_link
WHERE FALSE;

ROLLBACK;
