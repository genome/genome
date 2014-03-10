-- Deploy subject_pairing_permission
-- requires: subject_pairing

BEGIN;

REVOKE ALL ON TABLE subject.pairing FROM PUBLIC;
REVOKE ALL ON TABLE subject.pairing FROM genome;
GRANT ALL ON TABLE subject.pairing TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE subject.pairing TO "gms-user";

COMMIT;
