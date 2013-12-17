-- Revert subject_subject_permission

BEGIN;

REVOKE ALL ON TABLE subject.subject FROM "gms-user";

COMMIT;
