-- Revert subject_project_part_permission

BEGIN;

REVOKE ALL ON TABLE subject.project_part FROM "gms-user";

COMMIT;
