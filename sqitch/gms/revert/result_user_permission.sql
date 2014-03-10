-- Revert result_user_permission

BEGIN;

REVOKE ALL ON TABLE result."user" FROM "gms-user";

COMMIT;
