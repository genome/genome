-- Revert subject_misc_attribute_permission

BEGIN;

REVOKE ALL ON TABLE subject.misc_attribute FROM "gms-user";

COMMIT;
