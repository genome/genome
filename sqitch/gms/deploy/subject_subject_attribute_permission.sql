-- Deploy subject_subject_attribute_permission
-- requires: subject_subject_attribute

BEGIN;

REVOKE ALL ON TABLE subject.subject_attribute FROM PUBLIC;
REVOKE ALL ON TABLE subject.subject_attribute FROM genome;
GRANT ALL ON TABLE subject.subject_attribute TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE subject.subject_attribute TO "gms-user";

COMMIT;
