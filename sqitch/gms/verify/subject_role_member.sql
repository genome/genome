-- Verify subject_role_member

BEGIN;

SELECT user_email, role_id
FROM subject.role_member
WHERE FALSE;

ROLLBACK;
