-- Deploy subject_misc_note_permission
-- requires: subject_misc_note

BEGIN;

REVOKE ALL ON TABLE subject.misc_note FROM PUBLIC;
REVOKE ALL ON TABLE subject.misc_note FROM genome;
GRANT ALL ON TABLE subject.misc_note TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE subject.misc_note TO "gms-user";

COMMIT;
