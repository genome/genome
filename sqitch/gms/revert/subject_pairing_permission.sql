-- Revert subject_pairing_permission

BEGIN;

REVOKE ALL ON TABLE subject.pairing FROM "gms-user";

COMMIT;
