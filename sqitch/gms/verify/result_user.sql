-- Verify result_user

BEGIN;

SELECT id, software_result_id, user_id, user_class_name, label, active
FROM result."user"
WHERE FALSE;

ROLLBACK;
