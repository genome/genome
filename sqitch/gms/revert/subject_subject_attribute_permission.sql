-- Revert subject_subject_attribute_permission

BEGIN;

REVOKE ALL ON TABLE subject.subject_attribute FROM "gms-user";

COMMIT;
