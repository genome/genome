-- Revert subject_user_permission

BEGIN;

REVOKE ALL ON TABLE subject."user" FROM "gms-user";

COMMIT;
