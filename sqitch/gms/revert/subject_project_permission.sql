-- Revert subject_project_permission

BEGIN;

REVOKE ALL ON TABLE subject.project FROM "gms-user";

COMMIT;
