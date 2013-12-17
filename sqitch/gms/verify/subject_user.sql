-- Verify subject_user

BEGIN;

SELECT name, email, username
FROM subject."user"
WHERE FALSE;

ROLLBACK;
