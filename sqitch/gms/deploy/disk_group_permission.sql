-- Deploy disk_group_permission
-- requires: disk_group

BEGIN;

REVOKE ALL ON TABLE disk."group" FROM PUBLIC;
REVOKE ALL ON TABLE disk."group" FROM genome;
GRANT ALL ON TABLE disk."group" TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE disk."group" TO "gms-user";

COMMIT;
