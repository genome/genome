-- Revert subject_misc_note_permission

BEGIN;

REVOKE ALL ON TABLE subject.misc_note FROM "gms-user";

COMMIT;
