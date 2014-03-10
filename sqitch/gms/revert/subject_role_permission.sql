-- Revert subject_role_permission

BEGIN;

REVOKE ALL ON TABLE subject.role FROM "gms-user";

COMMIT;
