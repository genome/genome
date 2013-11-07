-- Revert subject_role_member_permission

BEGIN;

REVOKE ALL ON TABLE subject.role_member FROM "gms-user";

COMMIT;
