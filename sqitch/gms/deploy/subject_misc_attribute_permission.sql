-- Deploy subject_misc_attribute_permission
-- requires: subject_misc_attribute

BEGIN;

REVOKE ALL ON TABLE subject.misc_attribute FROM PUBLIC;
REVOKE ALL ON TABLE subject.misc_attribute FROM genome;
GRANT ALL ON TABLE subject.misc_attribute TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE subject.misc_attribute TO "gms-user";

COMMIT;
