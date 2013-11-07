-- Deploy subject_subject_permission
-- requires: subject_subject

BEGIN;

REVOKE ALL ON TABLE subject.subject FROM PUBLIC;
REVOKE ALL ON TABLE subject.subject FROM genome;
GRANT ALL ON TABLE subject.subject TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE subject.subject TO "gms-user";

COMMIT;
